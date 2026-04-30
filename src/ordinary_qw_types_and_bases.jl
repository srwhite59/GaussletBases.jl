"""
    OrdinaryCartesianOrbital3D

One orbital index in the paper-faithful Qiu-White residual-Gaussian reference
path.

The object records whether the orbital is a Cartesian gausslet product orbital
or a residual Gaussian together with the associated center and, for residual
Gaussians, the matched widths used in the MWG diagnostics. `owner_nucleus_index`
is `0` for unowned fixed/gausslet orbitals and the 1-based nucleus index for
atom-local residual Gaussians.
"""
struct OrdinaryCartesianOrbital3D
    index::Int
    kind::Symbol
    label::String
    x::Float64
    y::Float64
    z::Float64
    wx::Float64
    wy::Float64
    wz::Float64
    owner_nucleus_index::Int
end

function OrdinaryCartesianOrbital3D(
    index::Integer,
    kind::Symbol,
    label::AbstractString,
    x::Real,
    y::Real,
    z::Real,
    wx::Real,
    wy::Real,
    wz::Real,
)
    return OrdinaryCartesianOrbital3D(
        Int(index),
        kind,
        String(label),
        Float64(x),
        Float64(y),
        Float64(z),
        Float64(wx),
        Float64(wy),
        Float64(wz),
        0,
    )
end

const QiuWhiteHybridOrbital3D = OrdinaryCartesianOrbital3D

function Base.show(io::IO, orbital::OrdinaryCartesianOrbital3D)
    print(
        io,
        "OrdinaryCartesianOrbital3D(index=",
        orbital.index,
        ", kind=:",
        orbital.kind,
        ", label=\"",
        orbital.label,
        "\", center=(",
        orbital.x,
        ", ",
        orbital.y,
        ", ",
        orbital.z,
        ")",
    )
    if orbital.kind == :residual_gaussian
        print(
            io,
            ", widths=(",
            orbital.wx,
            ", ",
            orbital.wy,
            ", ",
            orbital.wz,
            ")",
        )
    end
    orbital.owner_nucleus_index == 0 || print(io, ", owner_nucleus_index=", orbital.owner_nucleus_index)
    print(io, ")")
end

"""
    OrdinaryCartesianOperators3D

Ordinary Cartesian operator payload for the current pure-Cartesian and
hybrid-residual construction routes.

This object keeps the final basis as the full 3D gausslet product basis plus
orthonormalized 3D residual Gaussians. The one-body matrices are built exactly
in the raw gausslet-plus-GTO space and transformed into the final basis when a
residual supplement is present, while the two-electron interaction stays in the
same two-index integral-diagonal approximation (IDA) representation used for
the gausslet channel.
"""
struct OrdinaryCartesianOperators3D{B,D}
    basis::B
    gaussian_data::D
    gausslet_backend::Symbol
    interaction_treatment::Symbol
    expansion::CoulombGaussianExpansion
    overlap::Matrix{Float64}
    one_body_hamiltonian::Matrix{Float64}
    interaction_matrix::Matrix{Float64}
    orbital_data::Vector{OrdinaryCartesianOrbital3D}
    gausslet_count::Int
    residual_count::Int
    raw_to_final::Matrix{Float64}
    residual_centers::Matrix{Float64}
    residual_widths::Matrix{Float64}
    residual_nucleus_indices::Vector{Int}
    nuclear_charges::Union{Nothing,Vector{Float64}}
    kinetic_one_body::Union{Nothing,Matrix{Float64}}
    nuclear_one_body_by_center::Union{Nothing,Vector{Matrix{Float64}}}
    nuclear_term_storage::Symbol
end

const QiuWhiteResidualGaussianOperators = OrdinaryCartesianOperators3D

function OrdinaryCartesianOperators3D(
    basis,
    gaussian_data,
    gausslet_backend,
    interaction_treatment,
    expansion,
    overlap,
    one_body_hamiltonian,
    interaction_matrix,
    orbital_data,
    gausslet_count,
    residual_count,
    raw_to_final,
    residual_centers,
    residual_widths,
)
    return OrdinaryCartesianOperators3D(
        basis,
        gaussian_data,
        gausslet_backend,
        interaction_treatment,
        expansion,
        overlap,
        one_body_hamiltonian,
        interaction_matrix,
        orbital_data,
        gausslet_count,
        residual_count,
        raw_to_final,
        residual_centers,
        residual_widths,
        zeros(Int, residual_count),
        nothing,
        nothing,
        nothing,
        :total_only,
    )
end

function OrdinaryCartesianOperators3D(
    basis,
    gaussian_data,
    gausslet_backend,
    interaction_treatment,
    expansion,
    overlap,
    one_body_hamiltonian,
    interaction_matrix,
    orbital_data,
    gausslet_count,
    residual_count,
    raw_to_final,
    residual_centers,
    residual_widths,
    nuclear_charges,
    kinetic_one_body,
    nuclear_one_body_by_center,
    nuclear_term_storage,
)
    return OrdinaryCartesianOperators3D(
        basis,
        gaussian_data,
        gausslet_backend,
        interaction_treatment,
        expansion,
        overlap,
        one_body_hamiltonian,
        interaction_matrix,
        orbital_data,
        gausslet_count,
        residual_count,
        raw_to_final,
        residual_centers,
        residual_widths,
        zeros(Int, residual_count),
        nuclear_charges,
        kinetic_one_body,
        nuclear_one_body_by_center,
        nuclear_term_storage,
    )
end

"""
    QWRGResidualSpaceDiagnostics

Focused overlap/rank diagnostics for the Qiu-White residual-space construction.

This records the raw supplement-space overlap spectrum, the residualized
supplement spectrum before the keep rule is applied, the diagnostic null-rank
threshold, the actual keep threshold used by the current route, and the final
kept/discarded counts.
"""
struct QWRGResidualSpaceDiagnostics
    gausslet_count::Int
    gaussian_count::Int
    raw_overlap_dimension::Int
    overlap_ga_size::Tuple{Int,Int}
    overlap_aa_size::Tuple{Int,Int}
    gausslet_overlap_error::Float64
    supplement_null_rank_tol::Float64
    supplement_numerical_rank::Int
    residual_null_rank_tol::Float64
    residual_numerical_rank::Int
    keep_policy::Symbol
    keep_tol::Float64
    accept_tol::Float64
    kept_count::Int
    discarded_count::Int
    kept_block_stabilization_null_tol::Float64
    kept_block_stabilization_correction_passes::Int
    kept_block_stabilization_clipped_count::Int
    kept_block_stabilization_dropped_count::Int
    kept_block_pre_stabilization_overlap_error::Float64
    kept_block_post_stabilization_overlap_error::Float64
    kept_block_pre_stabilization_symmetry_defect::Float64
    kept_block_post_stabilization_symmetry_defect::Float64
    kept_block_pre_stabilization_min_eigenvalue::Float64
    kept_block_pre_stabilization_max_eigenvalue::Float64
    kept_block_post_stabilization_min_eigenvalue::Float64
    kept_block_post_stabilization_max_eigenvalue::Float64
    kept_block_pre_stabilization_negative_count::Int
    kept_block_post_stabilization_negative_count::Int
    kept_block_pre_stabilization_near_null_count::Int
    kept_block_post_stabilization_near_null_count::Int
    supplement_overlap_eigenvalues::Vector{Float64}
    residual_overlap_eigenvalues::Vector{Float64}
    kept_indices::Vector{Int}
    discarded_indices::Vector{Int}
    kept_eigenvalues::Vector{Float64}
    discarded_eigenvalues::Vector{Float64}
end

function Base.show(io::IO, diagnostics::QWRGResidualSpaceDiagnostics)
    print(
        io,
        "QWRGResidualSpaceDiagnostics(ngausslet=",
        diagnostics.gausslet_count,
        ", ngaussian=",
        diagnostics.gaussian_count,
        ", residual_rank=",
        diagnostics.residual_numerical_rank,
        ", kept=",
        diagnostics.kept_count,
        ", keep_policy=:",
        diagnostics.keep_policy,
        ", keep_tol=",
        diagnostics.keep_tol,
        ", accept_tol=",
        diagnostics.accept_tol,
        ", stab_passes=",
        diagnostics.kept_block_stabilization_correction_passes,
        ", pre_stab_error=",
        diagnostics.kept_block_pre_stabilization_overlap_error,
        ", post_stab_error=",
        diagnostics.kept_block_post_stabilization_overlap_error,
        ")",
    )
end

"""
    BondAlignedDiatomicQWBasis3D

Narrow mixed-axis basis container for the first bond-aligned diatomic QW
reference route.

The first supported geometry family is a bond-aligned homonuclear diatomic:

- one combined multi-center mapping on the distinguished bond axis
- one shared single-center mapping on the two transverse axes
- one rectangular 3D product basis built from those three one-dimensional
  mapped bases
"""
abstract type AbstractBondAlignedOrdinaryQWBasis3D end

struct BondAlignedDiatomicQWBasis3D{B<:MappedUniformBasis} <: AbstractBondAlignedOrdinaryQWBasis3D
    bond_axis::Symbol
    basis_x::B
    basis_y::B
    basis_z::B
    nuclei::Vector{NTuple{3,Float64}}
    nuclear_charges::Vector{Float64}
    target_core_spacing::Float64
end

function Base.show(io::IO, basis::BondAlignedDiatomicQWBasis3D)
    print(
        io,
        "BondAlignedDiatomicQWBasis3D(bond_axis=:",
        basis.bond_axis,
        ", nx=",
        length(basis.basis_x),
        ", ny=",
        length(basis.basis_y),
        ", nz=",
        length(basis.basis_z),
        ", nuclei=",
        basis.nuclei,
        ", nuclear_charges=",
        basis.nuclear_charges,
        ", target_core_spacing=",
        basis.target_core_spacing,
        ")",
    )
end

"""
    BondAlignedHomonuclearChainQWBasis3D

Experimental ordinary QW basis container for a homonuclear linear chain on one
distinguished Cartesian axis.
"""
struct BondAlignedHomonuclearChainQWBasis3D{B<:MappedUniformBasis} <: AbstractBondAlignedOrdinaryQWBasis3D
    chain_axis::Symbol
    basis_x::B
    basis_y::B
    basis_z::B
    nuclei::Vector{NTuple{3,Float64}}
    nuclear_charges::Vector{Float64}
    chain_coordinates::Vector{Float64}
    target_core_spacing::Float64
end

function Base.show(io::IO, basis::BondAlignedHomonuclearChainQWBasis3D)
    print(
        io,
        "BondAlignedHomonuclearChainQWBasis3D(chain_axis=:",
        basis.chain_axis,
        ", natoms=",
        length(basis.nuclei),
        ", nx=",
        length(basis.basis_x),
        ", ny=",
        length(basis.basis_y),
        ", nz=",
        length(basis.basis_z),
        ", nuclear_charges=",
        basis.nuclear_charges,
        ", target_core_spacing=",
        basis.target_core_spacing,
        ")",
    )
end

"""
    AxisAlignedHomonuclearSquareLatticeQWBasis3D

Experimental ordinary QW basis container for an axis-aligned homonuclear
`n × n` square lattice in the `xy` plane.
"""
struct AxisAlignedHomonuclearSquareLatticeQWBasis3D{B<:MappedUniformBasis} <: AbstractBondAlignedOrdinaryQWBasis3D
    lattice_size::Int
    basis_x::B
    basis_y::B
    basis_z::B
    nuclei::Vector{NTuple{3,Float64}}
    nuclear_charges::Vector{Float64}
    x_coordinates::Vector{Float64}
    y_coordinates::Vector{Float64}
    target_core_spacing::Float64
end

function Base.show(io::IO, basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    print(
        io,
        "AxisAlignedHomonuclearSquareLatticeQWBasis3D(n=",
        basis.lattice_size,
        ", nx=",
        length(basis.basis_x),
        ", ny=",
        length(basis.basis_y),
        ", nz=",
        length(basis.basis_z),
        ", nuclear_charges=",
        basis.nuclear_charges,
        ", target_core_spacing=",
        basis.target_core_spacing,
        ")",
    )
end

function _qwrg_axis_coordinate(
    point::NTuple{3,Float64},
    axis::Symbol,
)
    axis == :x && return point[1]
    axis == :y && return point[2]
    axis == :z && return point[3]
    throw(ArgumentError("axis must be :x, :y, or :z"))
end

function _qwrg_mapped_odd_count_for_extent(
    mapping_value::AbstractCoordinateMapping,
    xmax::Real;
    reference_spacing::Real = 1.0,
)
    xmax_value = Float64(xmax)
    spacing_value = Float64(reference_spacing)
    xmax_value > 0.0 || throw(ArgumentError("mapped extent helper requires xmax > 0"))
    spacing_value > 0.0 || throw(ArgumentError("mapped extent helper requires reference_spacing > 0"))
    uedge = uofx(mapping_value, xmax_value)
    count = 2 * ceil(Int, uedge / spacing_value) + 1
    isodd(count) || throw(ArgumentError("mapped extent helper must produce an odd count"))
    return count
end

function _qwrg_mapped_odd_count_for_interval(
    mapping_value::AbstractCoordinateMapping,
    xmin::Real,
    xmax::Real;
    reference_spacing::Real = 1.0,
)
    xmin_value = Float64(xmin)
    xmax_value = Float64(xmax)
    xmin_value < xmax_value || throw(ArgumentError("mapped interval helper requires xmin < xmax"))
    spacing_value = Float64(reference_spacing)
    spacing_value > 0.0 || throw(ArgumentError("mapped interval helper requires reference_spacing > 0"))
    uextent = max(abs(uofx(mapping_value, xmin_value)), abs(uofx(mapping_value, xmax_value)))
    count = 2 * ceil(Int, uextent / spacing_value) + 1
    isodd(count) || throw(ArgumentError("mapped interval helper must produce an odd count"))
    return count
end

function _qwrg_centered_chain_coordinates(
    natoms::Integer,
    spacing::Real,
)
    natoms > 0 || throw(ArgumentError("homonuclear chain helper requires natoms > 0"))
    spacing_value = Float64(spacing)
    spacing_value > 0.0 || throw(ArgumentError("homonuclear chain helper requires spacing > 0"))
    midpoint = 0.5 * (Int(natoms) + 1)
    return Float64[(index - midpoint) * spacing_value for index in 1:Int(natoms)]
end

function _qwrg_validate_chain_coordinates(chain_coordinates::AbstractVector{<:Real})
    length(chain_coordinates) > 0 || throw(ArgumentError("homonuclear chain basis requires at least one nuclear coordinate"))
    values = Float64[Float64(value) for value in chain_coordinates]
    issorted(values) || throw(ArgumentError("explicit homonuclear chain coordinates must be ordered along the chain axis"))
    for index in 2:length(values)
        values[index] > values[index - 1] || throw(ArgumentError("explicit homonuclear chain coordinates must be strictly increasing"))
    end
    return values
end

function _qwrg_validate_square_lattice_axis_coordinates(
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    atol::Float64 = 1.0e-10,
    rtol::Float64 = 1.0e-8,
)
    x_values = _qwrg_validate_chain_coordinates(x_coordinates)
    y_values = _qwrg_validate_chain_coordinates(y_coordinates)
    length(x_values) == length(y_values) || throw(
        ArgumentError("square lattice explicit coordinates require equally sized x and y coordinate lists"),
    )
    length(x_values) > 0 || throw(ArgumentError("square lattice explicit coordinates require at least one site per axis"))

    if length(x_values) > 1
        x_steps = diff(x_values)
        y_steps = diff(y_values)
        x_step = first(x_steps)
        y_step = first(y_steps)
        all(step -> isapprox(step, x_step; atol = atol, rtol = rtol), x_steps) || throw(
            ArgumentError("square lattice explicit x coordinates must be uniformly spaced"),
        )
        all(step -> isapprox(step, y_step; atol = atol, rtol = rtol), y_steps) || throw(
            ArgumentError("square lattice explicit y coordinates must be uniformly spaced"),
        )
        isapprox(x_step, y_step; atol = atol, rtol = rtol) || throw(
            ArgumentError("square lattice explicit x and y coordinates must use the same lattice spacing"),
        )
    end

    return x_values, y_values
end

function _qwrg_homonuclear_chain_nuclei(
    chain_coordinates::AbstractVector{<:Real},
    chain_axis::Symbol,
)
    coordinates = _qwrg_validate_chain_coordinates(chain_coordinates)
    chain_axis == :x && return NTuple{3,Float64}[(value, 0.0, 0.0) for value in coordinates]
    chain_axis == :y && return NTuple{3,Float64}[(0.0, value, 0.0) for value in coordinates]
    chain_axis == :z && return NTuple{3,Float64}[(0.0, 0.0, value) for value in coordinates]
    throw(ArgumentError("chain_axis must be :x, :y, or :z"))
end

function _qwrg_homonuclear_square_lattice_nuclei(
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real},
)
    x_values, y_values = _qwrg_validate_square_lattice_axis_coordinates(
        x_coordinates,
        y_coordinates,
    )
    return NTuple{3,Float64}[(x, y, 0.0) for y in y_values for x in x_values]
end

function _qwrg_bond_aligned_homonuclear_nuclei(
    bond_length::Real,
    bond_axis::Symbol,
)
    half = 0.5 * Float64(bond_length)
    bond_axis == :x && return [(-half, 0.0, 0.0), (half, 0.0, 0.0)]
    bond_axis == :y && return [(0.0, -half, 0.0), (0.0, half, 0.0)]
    bond_axis == :z && return [(0.0, 0.0, -half), (0.0, 0.0, half)]
    throw(ArgumentError("bond_axis must be :x, :y, or :z"))
end

function _qwrg_bond_aligned_diatomic_nuclei(
    bond_length::Real,
    bond_axis::Symbol,
)
    return _qwrg_bond_aligned_homonuclear_nuclei(bond_length, bond_axis)
end

function _qwrg_tighter_heavier_index(
    core_spacings::NTuple{2,Float64},
    nuclear_charges::NTuple{2,Float64},
)
    if core_spacings[1] < core_spacings[2]
        return 1
    elseif core_spacings[2] < core_spacings[1]
        return 2
    elseif nuclear_charges[1] >= nuclear_charges[2]
        return 1
    else
        return 2
    end
end

"""
    bond_aligned_homonuclear_qw_basis(; ...)

Build the first bond-aligned homonuclear diatomic 3D product basis for the
ordinary QW reference line.

The bond axis uses a combined multi-center inverse-sqrt-density mapping, while
the two transverse axes share a single-center inverse-sqrt mapping at the
common transverse projection.
"""
function bond_aligned_homonuclear_qw_basis(;
    family = :G10,
    bond_length::Real,
    core_spacing::Real = 0.5,
    xmax_parallel::Real = 8.0,
    xmax_transverse::Real = 6.0,
    bond_axis::Symbol = :z,
    nuclear_charge::Real = 1.0,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
)
    core_spacing_value = Float64(core_spacing)
    nuclear_charge_value = Float64(nuclear_charge)
    core_spacing_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_qw_basis requires core_spacing > 0"))
    nuclear_charge_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_qw_basis requires nuclear_charge > 0"))

    nuclei = _qwrg_bond_aligned_homonuclear_nuclei(bond_length, bond_axis)
    parallel_centers = Float64[_qwrg_axis_coordinate(nucleus, bond_axis) for nucleus in nuclei]
    core_range = sqrt(core_spacing_value / nuclear_charge_value)

    # Alg Nested-Diatomic-Map step 4 and 6: use a combined inverse-sqrt map on
    # the bond axis and one shared single-center map on the transverse axes.
    # See docs/src/algorithms/cartesian_nested_diatomic_coordinate_distortion.md.
    parallel_mapping = fit_combined_invsqrt_mapping(
        centers = parallel_centers,
        core_ranges = fill(core_range, length(parallel_centers)),
        target_spacings = fill(core_spacing_value, length(parallel_centers)),
        tail_spacing = tail_spacing,
    )
    transverse_mapping = fit_combined_invsqrt_mapping(
        centers = [0.0],
        core_ranges = [core_range],
        target_spacings = [core_spacing_value],
        tail_spacing = tail_spacing,
    )

    count_parallel = _qwrg_mapped_odd_count_for_extent(
        parallel_mapping,
        xmax_parallel;
        reference_spacing = reference_spacing,
    )
    count_transverse = _qwrg_mapped_odd_count_for_extent(
        transverse_mapping,
        xmax_transverse;
        reference_spacing = reference_spacing,
    )

    parallel_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_parallel,
        mapping = parallel_mapping,
        reference_spacing = reference_spacing,
    ))
    transverse_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_transverse,
        mapping = transverse_mapping,
        reference_spacing = reference_spacing,
    ))

    basis_x = bond_axis == :x ? parallel_basis : transverse_basis
    basis_y = bond_axis == :y ? parallel_basis : transverse_basis
    basis_z = bond_axis == :z ? parallel_basis : transverse_basis
    return BondAlignedDiatomicQWBasis3D(
        bond_axis,
        basis_x,
        basis_y,
        basis_z,
        nuclei,
        fill(nuclear_charge_value, length(nuclei)),
        core_spacing_value,
    )
end

"""
    bond_aligned_homonuclear_chain_qw_basis(; ...)

Build the first experimental ordinary QW basis for a homonuclear linear chain.

The chain axis uses a combined inverse-sqrt map fitted to all nuclear centers
on that axis. The transverse axes share one single-center combined inverse-sqrt
map at the common transverse projection.
"""
function bond_aligned_homonuclear_chain_qw_basis(;
    family = :G10,
    natoms::Union{Nothing,Integer} = nothing,
    spacing::Union{Nothing,Real} = nothing,
    chain_coordinates::Union{Nothing,AbstractVector{<:Real}} = nothing,
    core_spacing::Real = 0.5,
    xmax_parallel::Real = 4.0,
    xmax_transverse::Real = 4.0,
    chain_axis::Symbol = :z,
    nuclear_charge::Real = 1.0,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
)
    helper_mode = natoms !== nothing || spacing !== nothing
    explicit_mode = chain_coordinates !== nothing
    helper_mode ⊻ explicit_mode || throw(
        ArgumentError("bond_aligned_homonuclear_chain_qw_basis requires either natoms+spacing or explicit chain_coordinates"),
    )
    helper_mode && (natoms === nothing || spacing === nothing) && throw(
        ArgumentError("bond_aligned_homonuclear_chain_qw_basis requires both natoms and spacing in helper mode"),
    )

    core_spacing_value = Float64(core_spacing)
    nuclear_charge_value = Float64(nuclear_charge)
    xmax_parallel_value = Float64(xmax_parallel)
    xmax_transverse_value = Float64(xmax_transverse)
    core_spacing_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_chain_qw_basis requires core_spacing > 0"))
    nuclear_charge_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_chain_qw_basis requires nuclear_charge > 0"))
    xmax_parallel_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_chain_qw_basis requires xmax_parallel > 0"))
    xmax_transverse_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_chain_qw_basis requires xmax_transverse > 0"))

    chain_coordinate_values =
        explicit_mode ? _qwrg_validate_chain_coordinates(chain_coordinates) :
        _qwrg_centered_chain_coordinates(natoms, spacing)
    nuclei = _qwrg_homonuclear_chain_nuclei(chain_coordinate_values, chain_axis)
    core_range = sqrt(core_spacing_value / nuclear_charge_value)

    parallel_mapping = fit_combined_invsqrt_mapping(
        centers = chain_coordinate_values,
        core_ranges = fill(core_range, length(chain_coordinate_values)),
        target_spacings = fill(core_spacing_value, length(chain_coordinate_values)),
        tail_spacing = tail_spacing,
    )
    transverse_mapping = fit_combined_invsqrt_mapping(
        centers = [0.0],
        core_ranges = [core_range],
        target_spacings = [core_spacing_value],
        tail_spacing = tail_spacing,
    )

    count_parallel = _qwrg_mapped_odd_count_for_interval(
        parallel_mapping,
        minimum(chain_coordinate_values) - xmax_parallel_value,
        maximum(chain_coordinate_values) + xmax_parallel_value;
        reference_spacing = reference_spacing,
    )
    count_transverse = _qwrg_mapped_odd_count_for_extent(
        transverse_mapping,
        xmax_transverse_value;
        reference_spacing = reference_spacing,
    )

    parallel_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_parallel,
        mapping = parallel_mapping,
        reference_spacing = reference_spacing,
    ))
    transverse_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_transverse,
        mapping = transverse_mapping,
        reference_spacing = reference_spacing,
    ))

    basis_x = chain_axis == :x ? parallel_basis : transverse_basis
    basis_y = chain_axis == :y ? parallel_basis : transverse_basis
    basis_z = chain_axis == :z ? parallel_basis : transverse_basis
    return BondAlignedHomonuclearChainQWBasis3D(
        chain_axis,
        basis_x,
        basis_y,
        basis_z,
        nuclei,
        fill(nuclear_charge_value, length(nuclei)),
        chain_coordinate_values,
        core_spacing_value,
    )
end

"""
    axis_aligned_homonuclear_square_lattice_qw_basis(; ...)

Build the first experimental ordinary QW basis for an axis-aligned homonuclear
`n × n` square lattice in the `xy` plane.

The in-plane `x` and `y` axes each use a combined inverse-sqrt map fitted to
the unique lattice coordinates on that axis. The transverse `z` axis uses one
shared single-center combined inverse-sqrt map at the common lattice plane.
"""
function axis_aligned_homonuclear_square_lattice_qw_basis(;
    family = :G10,
    n::Union{Nothing,Integer} = nothing,
    spacing::Union{Nothing,Real} = nothing,
    x_coordinates::Union{Nothing,AbstractVector{<:Real}} = nothing,
    y_coordinates::Union{Nothing,AbstractVector{<:Real}} = nothing,
    core_spacing::Real = 0.5,
    xmax_in_plane::Real = 4.0,
    xmax_transverse::Real = 4.0,
    nuclear_charge::Real = 1.0,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
)
    helper_mode = n !== nothing || spacing !== nothing
    explicit_mode = x_coordinates !== nothing || y_coordinates !== nothing
    helper_mode ⊻ explicit_mode || throw(
        ArgumentError("axis_aligned_homonuclear_square_lattice_qw_basis requires either n+spacing or explicit x_coordinates+y_coordinates"),
    )
    helper_mode && (n === nothing || spacing === nothing) && throw(
        ArgumentError("axis_aligned_homonuclear_square_lattice_qw_basis requires both n and spacing in helper mode"),
    )
    explicit_mode && (x_coordinates === nothing || y_coordinates === nothing) && throw(
        ArgumentError("axis_aligned_homonuclear_square_lattice_qw_basis requires both x_coordinates and y_coordinates in explicit mode"),
    )

    core_spacing_value = Float64(core_spacing)
    nuclear_charge_value = Float64(nuclear_charge)
    xmax_in_plane_value = Float64(xmax_in_plane)
    xmax_transverse_value = Float64(xmax_transverse)
    core_spacing_value > 0.0 || throw(ArgumentError("axis_aligned_homonuclear_square_lattice_qw_basis requires core_spacing > 0"))
    nuclear_charge_value > 0.0 || throw(ArgumentError("axis_aligned_homonuclear_square_lattice_qw_basis requires nuclear_charge > 0"))
    xmax_in_plane_value > 0.0 || throw(ArgumentError("axis_aligned_homonuclear_square_lattice_qw_basis requires xmax_in_plane > 0"))
    xmax_transverse_value > 0.0 || throw(ArgumentError("axis_aligned_homonuclear_square_lattice_qw_basis requires xmax_transverse > 0"))

    x_coordinate_values, y_coordinate_values =
        explicit_mode ? _qwrg_validate_square_lattice_axis_coordinates(
            x_coordinates,
            y_coordinates,
        ) :
        (
            _qwrg_centered_chain_coordinates(n, spacing),
            _qwrg_centered_chain_coordinates(n, spacing),
        )
    lattice_size = length(x_coordinate_values)
    nuclei = _qwrg_homonuclear_square_lattice_nuclei(x_coordinate_values, y_coordinate_values)
    core_range = sqrt(core_spacing_value / nuclear_charge_value)

    x_mapping = fit_combined_invsqrt_mapping(
        centers = x_coordinate_values,
        core_ranges = fill(core_range, lattice_size),
        target_spacings = fill(core_spacing_value, lattice_size),
        tail_spacing = tail_spacing,
    )
    y_mapping = fit_combined_invsqrt_mapping(
        centers = y_coordinate_values,
        core_ranges = fill(core_range, lattice_size),
        target_spacings = fill(core_spacing_value, lattice_size),
        tail_spacing = tail_spacing,
    )
    z_mapping = fit_combined_invsqrt_mapping(
        centers = [0.0],
        core_ranges = [core_range],
        target_spacings = [core_spacing_value],
        tail_spacing = tail_spacing,
    )

    count_x = _qwrg_mapped_odd_count_for_interval(
        x_mapping,
        minimum(x_coordinate_values) - xmax_in_plane_value,
        maximum(x_coordinate_values) + xmax_in_plane_value;
        reference_spacing = reference_spacing,
    )
    count_y = _qwrg_mapped_odd_count_for_interval(
        y_mapping,
        minimum(y_coordinate_values) - xmax_in_plane_value,
        maximum(y_coordinate_values) + xmax_in_plane_value;
        reference_spacing = reference_spacing,
    )
    count_z = _qwrg_mapped_odd_count_for_extent(
        z_mapping,
        xmax_transverse_value;
        reference_spacing = reference_spacing,
    )

    basis_x = build_basis(MappedUniformBasisSpec(
        family;
        count = count_x,
        mapping = x_mapping,
        reference_spacing = reference_spacing,
    ))
    basis_y = build_basis(MappedUniformBasisSpec(
        family;
        count = count_y,
        mapping = y_mapping,
        reference_spacing = reference_spacing,
    ))
    basis_z = build_basis(MappedUniformBasisSpec(
        family;
        count = count_z,
        mapping = z_mapping,
        reference_spacing = reference_spacing,
    ))

    return AxisAlignedHomonuclearSquareLatticeQWBasis3D(
        lattice_size,
        basis_x,
        basis_y,
        basis_z,
        nuclei,
        fill(nuclear_charge_value, length(nuclei)),
        x_coordinate_values,
        y_coordinate_values,
        core_spacing_value,
    )
end

"""
    bond_aligned_heteronuclear_qw_basis(; ...)

Build the first bond-aligned heteronuclear diatomic 3D product basis for the
ordinary QW reference line.

The bond axis stays on the current combined inverse-sqrt-density family, but
now with explicit per-atom local spacing targets. The transverse axes still use
one shared single-center inverse-sqrt mapping at the common transverse
projection, controlled by the tighter/heavier side.
"""
function bond_aligned_heteronuclear_qw_basis(;
    family = :G10,
    atoms::Tuple{<:AbstractString,<:AbstractString},
    bond_length::Real,
    core_spacings::Tuple{<:Real,<:Real},
    nuclear_charges::Tuple{<:Real,<:Real},
    xmax_parallel::Real = 8.0,
    xmax_transverse::Real = 6.0,
    transverse_core_spacing::Union{Nothing,Real} = nothing,
    bond_axis::Symbol = :z,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
)
    atom_a = String(atoms[1])
    atom_b = String(atoms[2])
    !isempty(atom_a) || throw(ArgumentError("bond_aligned_heteronuclear_qw_basis requires a nonempty first atom label"))
    !isempty(atom_b) || throw(ArgumentError("bond_aligned_heteronuclear_qw_basis requires a nonempty second atom label"))
    core_spacing_values = (Float64(core_spacings[1]), Float64(core_spacings[2]))
    nuclear_charge_values = (Float64(nuclear_charges[1]), Float64(nuclear_charges[2]))
    all(value -> value > 0.0, core_spacing_values) || throw(
        ArgumentError("bond_aligned_heteronuclear_qw_basis requires positive per-atom core spacings"),
    )
    all(value -> value > 0.0, nuclear_charge_values) || throw(
        ArgumentError("bond_aligned_heteronuclear_qw_basis requires positive nuclear charges"),
    )

    nuclei = _qwrg_bond_aligned_diatomic_nuclei(bond_length, bond_axis)
    parallel_centers = Float64[_qwrg_axis_coordinate(nucleus, bond_axis) for nucleus in nuclei]
    parallel_core_ranges = [
        sqrt(core_spacing_values[index] / nuclear_charge_values[index]) for index in 1:2
    ]
    transverse_index = _qwrg_tighter_heavier_index(core_spacing_values, nuclear_charge_values)
    transverse_spacing_value =
        transverse_core_spacing === nothing ? core_spacing_values[transverse_index] :
        Float64(transverse_core_spacing)
    transverse_spacing_value > 0.0 || throw(
        ArgumentError("bond_aligned_heteronuclear_qw_basis requires transverse_core_spacing > 0 when provided"),
    )
    transverse_core_range = sqrt(
        transverse_spacing_value / nuclear_charge_values[transverse_index],
    )

    # Alg Nested-Diatomic-Map step 4 and 6: keep the first heteronuclear bond
    # axis on the same combined inverse-sqrt family with explicit per-atom
    # spacing targets, and let the tighter/heavier side control the shared
    # transverse map.
    # See docs/src/algorithms/cartesian_nested_diatomic_coordinate_distortion.md.
    parallel_mapping = fit_combined_invsqrt_mapping(
        centers = parallel_centers,
        core_ranges = parallel_core_ranges,
        target_spacings = collect(core_spacing_values),
        tail_spacing = tail_spacing,
    )
    transverse_mapping = fit_combined_invsqrt_mapping(
        centers = [0.0],
        core_ranges = [transverse_core_range],
        target_spacings = [transverse_spacing_value],
        tail_spacing = tail_spacing,
    )

    count_parallel = _qwrg_mapped_odd_count_for_extent(
        parallel_mapping,
        xmax_parallel;
        reference_spacing = reference_spacing,
    )
    count_transverse = _qwrg_mapped_odd_count_for_extent(
        transverse_mapping,
        xmax_transverse;
        reference_spacing = reference_spacing,
    )

    parallel_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_parallel,
        mapping = parallel_mapping,
        reference_spacing = reference_spacing,
    ))
    transverse_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_transverse,
        mapping = transverse_mapping,
        reference_spacing = reference_spacing,
    ))

    basis_x = bond_axis == :x ? parallel_basis : transverse_basis
    basis_y = bond_axis == :y ? parallel_basis : transverse_basis
    basis_z = bond_axis == :z ? parallel_basis : transverse_basis
    return BondAlignedDiatomicQWBasis3D(
        bond_axis,
        basis_x,
        basis_y,
        basis_z,
        nuclei,
        collect(nuclear_charge_values),
        transverse_spacing_value,
    )
end

function _qwrg_basis_for_axis(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
    axis::Symbol,
)
    axis == :x && return basis.basis_x
    axis == :y && return basis.basis_y
    axis == :z && return basis.basis_z
    throw(ArgumentError("axis must be :x, :y, or :z"))
end

function _qwrg_bond_axis_local_spacings(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
)
    axis_symbol = basis isa BondAlignedDiatomicQWBasis3D ? basis.bond_axis : basis.chain_axis
    axis_basis = _qwrg_basis_for_axis(basis, axis_symbol)
    axis_mapping = mapping(axis_basis)
    return Float64[
        1.0 / dudx(axis_mapping, _qwrg_axis_coordinate(nucleus, axis_symbol)) for
        nucleus in basis.nuclei
    ]
end

function _qwrg_bond_axis_order(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
)
    axis_symbol = basis isa BondAlignedDiatomicQWBasis3D ? basis.bond_axis : basis.chain_axis
    coordinates = Float64[_qwrg_axis_coordinate(nucleus, axis_symbol) for nucleus in basis.nuclei]
    return sortperm(coordinates), coordinates
end

function bond_aligned_homonuclear_chain_geometry_diagnostics(
    basis::BondAlignedHomonuclearChainQWBasis3D,
)
    axis_basis = _qwrg_basis_for_axis(basis, basis.chain_axis)
    transverse_axis = basis.chain_axis == :x ? :y : :x
    transverse_basis = _qwrg_basis_for_axis(basis, transverse_axis)
    chain_coordinates = copy(basis.chain_coordinates)
    midpoint_spacings = if length(chain_coordinates) <= 1
        Float64[]
    else
        axis_mapping = mapping(axis_basis)
        Float64[
            1.0 / dudx(axis_mapping, 0.5 * (chain_coordinates[index] + chain_coordinates[index + 1])) for
            index in 1:(length(chain_coordinates) - 1)
        ]
    end
    return (
        chain_axis = basis.chain_axis,
        natoms = length(basis.nuclei),
        chain_coordinates = chain_coordinates,
        chain_axis_centers = centers(axis_basis),
        transverse_centers = centers(transverse_basis),
        local_spacings_at_nuclei = _qwrg_bond_axis_local_spacings(basis),
        local_spacings_at_midpoints = midpoint_spacings,
        axis_mapping_kind = typeof(mapping(axis_basis)),
        transverse_mapping_kind = typeof(mapping(transverse_basis)),
        axis_monotone = all(diff(centers(axis_basis)) .> 0.0),
        axis_center_symmetry_error = maximum(abs.(centers(axis_basis) .+ reverse(centers(axis_basis)))),
    )
end

function axis_aligned_homonuclear_square_lattice_geometry_diagnostics(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D,
)
    x_mapping = mapping(basis.basis_x)
    y_mapping = mapping(basis.basis_y)
    z_mapping = mapping(basis.basis_z)
    x_centers = centers(basis.basis_x)
    y_centers = centers(basis.basis_y)
    z_centers = centers(basis.basis_z)
    x_midpoints = if length(basis.x_coordinates) <= 1
        Float64[]
    else
        Float64[
            0.5 * (basis.x_coordinates[index] + basis.x_coordinates[index + 1]) for
            index in 1:(length(basis.x_coordinates) - 1)
        ]
    end
    y_midpoints = if length(basis.y_coordinates) <= 1
        Float64[]
    else
        Float64[
            0.5 * (basis.y_coordinates[index] + basis.y_coordinates[index + 1]) for
            index in 1:(length(basis.y_coordinates) - 1)
        ]
    end

    return (
        lattice_size = basis.lattice_size,
        natoms = length(basis.nuclei),
        x_coordinates = copy(basis.x_coordinates),
        y_coordinates = copy(basis.y_coordinates),
        x_axis_centers = x_centers,
        y_axis_centers = y_centers,
        z_axis_centers = z_centers,
        local_spacings_at_x_coordinates = Float64[
            1.0 / dudx(x_mapping, coordinate) for coordinate in basis.x_coordinates
        ],
        local_spacings_at_y_coordinates = Float64[
            1.0 / dudx(y_mapping, coordinate) for coordinate in basis.y_coordinates
        ],
        representative_midpoint_spacings_x = Float64[
            1.0 / dudx(x_mapping, coordinate) for coordinate in x_midpoints
        ],
        representative_midpoint_spacings_y = Float64[
            1.0 / dudx(y_mapping, coordinate) for coordinate in y_midpoints
        ],
        local_spacing_at_plane_center_x = 1.0 / dudx(x_mapping, 0.0),
        local_spacing_at_plane_center_y = 1.0 / dudx(y_mapping, 0.0),
        local_spacing_at_plane_center_z = 1.0 / dudx(z_mapping, 0.0),
        x_mapping_kind = typeof(x_mapping),
        y_mapping_kind = typeof(y_mapping),
        z_mapping_kind = typeof(z_mapping),
        x_axis_monotone = all(diff(x_centers) .> 0.0),
        y_axis_monotone = all(diff(y_centers) .> 0.0),
        z_axis_monotone = all(diff(z_centers) .> 0.0),
        x_axis_center_symmetry_error = maximum(abs.(x_centers .+ reverse(x_centers))),
        y_axis_center_symmetry_error = maximum(abs.(y_centers .+ reverse(y_centers))),
        xy_axis_center_match_error = length(x_centers) == length(y_centers) ?
            maximum(abs.(x_centers .- y_centers)) : Inf,
    )
end

function Base.show(io::IO, operators::OrdinaryCartesianOperators3D)
    print(
        io,
        "OrdinaryCartesianOperators3D(gausslet_backend=:",
        operators.gausslet_backend,
        ", interaction=:",
        operators.interaction_treatment,
        ", ngausslet=",
        operators.gausslet_count,
        ", nresidual=",
        operators.residual_count,
        ", nuclear_terms=:",
        operators.nuclear_term_storage,
        ", reference=true)",
    )
end

orbitals(operators::OrdinaryCartesianOperators3D) = operators.orbital_data
