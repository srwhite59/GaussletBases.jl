"""
    CartesianBasisMetadata3D

Minimal public metadata bundle identifying a three-dimensional Cartesian basis.

This is the Cartesian analogue of `BasisMetadata1D`. It records:

- the high-level basis kind, such as `:direct_product`, `:nested_fixed_block`,
  or `:hybrid_residual`
- the axis-level one-dimensional metadata reused from the existing 1D contract
- the defining parent-space contract and final basis size
- light basis-level geometric data such as labels, centers, and working-box
  metadata when that construction metadata already exists

The object is basis-defining metadata only. It is not an operator cache.
"""
struct CartesianBasisMetadata3D{AT <: NamedTuple, RT <: NamedTuple}
    basis_kind::Symbol
    axis_sharing::Symbol
    axis_metadata::AT
    parent_kind::Symbol
    parent_axis_counts::NTuple{3,Int}
    parent_dimension::Int
    final_dimension::Int
    working_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    basis_labels::Vector{String}
    basis_centers::Matrix{Float64}
    route_metadata::RT
end

function Base.show(io::IO, metadata::CartesianBasisMetadata3D)
    print(
        io,
        "CartesianBasisMetadata3D(kind=:",
        metadata.basis_kind,
        ", nfinal=",
        metadata.final_dimension,
        ", parent_kind=:",
        metadata.parent_kind,
        ", axis_sharing=:",
        metadata.axis_sharing,
        ")",
    )
end

"""
    CartesianBasisRepresentation3D

Public in-memory three-dimensional Cartesian basis representation.

The representation stores:

- `metadata`
- explicit 1D axis representations reused from the existing 1D contract
- the contraction from the stored parent/raw basis to the final basis
- parent-space labels and centers
- support-space indexing data when the construction is support-local
- optional parent-data sidecars needed to make mixed raw-space contracts
  explicit, for example the Cartesian parent representation and supplement
  orbital metadata in the QW residual-Gaussian route

This first-pass contract is representation-only. It deliberately avoids any
operator cache beyond the basis-defining contraction data itself.
"""
struct CartesianBasisRepresentation3D{
    MT <: CartesianBasisMetadata3D,
    AT <: NamedTuple,
    PT <: NamedTuple,
}
    metadata::MT
    axis_representations::AT
    contraction_kind::Symbol
    coefficient_matrix::Union{Nothing,Matrix{Float64}}
    parent_labels::Vector{String}
    parent_centers::Matrix{Float64}
    support_indices::Union{Nothing,Vector{Int}}
    support_states::Union{Nothing,Vector{NTuple{3,Int}}}
    parent_data::PT
end

function Base.show(io::IO, representation::CartesianBasisRepresentation3D)
    print(
        io,
        "CartesianBasisRepresentation3D(kind=:",
        representation.metadata.basis_kind,
        ", nfinal=",
        representation.metadata.final_dimension,
        ", parent_kind=:",
        representation.metadata.parent_kind,
        ", contraction=:",
        representation.contraction_kind,
        ")",
    )
end

basis_metadata(representation::CartesianBasisRepresentation3D) = representation.metadata

function _cartesian_product_states(axis_counts::NTuple{3,Int})
    nx, ny, nz = axis_counts
    states = Vector{NTuple{3,Int}}(undef, nx * ny * nz)
    index = 0
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        index += 1
        states[index] = (ix, iy, iz)
    end
    return states
end

function _cartesian_product_cross_overlap(
    left_states::AbstractVector{<:NTuple{3,Int}},
    right_states::AbstractVector{<:NTuple{3,Int}},
    overlap_x::AbstractMatrix{<:Real},
    overlap_y::AbstractMatrix{<:Real},
    overlap_z::AbstractMatrix{<:Real},
)
    matrix = zeros(Float64, length(left_states), length(right_states))
    for (row, (ix_left, iy_left, iz_left)) in pairs(left_states)
        for (column, (ix_right, iy_right, iz_right)) in pairs(right_states)
            matrix[row, column] =
                overlap_x[ix_left, ix_right] *
                overlap_y[iy_left, iy_right] *
                overlap_z[iz_left, iz_right]
        end
    end
    return matrix
end

function _cartesian_same_parent_raw_identity(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    left.metadata.parent_kind == :cartesian_product_basis || return false
    right.metadata.parent_kind == :cartesian_product_basis || return false
    left.metadata.parent_axis_counts == right.metadata.parent_axis_counts || return false
    isequal(left.metadata.axis_metadata, right.metadata.axis_metadata) || return false
    isequal(left.parent_labels, right.parent_labels) || return false
    isequal(left.parent_centers, right.parent_centers) || return false
    return true
end

function _cartesian_supported_exact_parent_overlap(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    return left.metadata.parent_kind == :cartesian_product_basis &&
           right.metadata.parent_kind == :cartesian_product_basis
end

function _cartesian_cross_overlap_error(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    if left.metadata.parent_kind == :cartesian_plus_supplement_raw ||
       right.metadata.parent_kind == :cartesian_plus_supplement_raw
        return ArgumentError(
            "exact Cartesian cross overlap is not yet implemented for hybrid/QW residual representations because the current public representation does not yet carry explicit exact raw cartesian-supplement cross-overlap metadata",
        )
    end
    return ArgumentError(
        "exact Cartesian cross overlap is not yet implemented for parent kinds :$(left.metadata.parent_kind) and :$(right.metadata.parent_kind)",
    )
end

function _cartesian_parent_state_basis(
    representation::CartesianBasisRepresentation3D,
)
    if representation.contraction_kind == :dense &&
       representation.support_indices !== nothing &&
       representation.support_states !== nothing
        return (
            states = representation.support_states,
            coefficients = Matrix{Float64}(representation.coefficient_matrix[representation.support_indices, :]),
        )
    elseif representation.contraction_kind == :dense
        return (
            states = _cartesian_product_states(representation.metadata.parent_axis_counts),
            coefficients = Matrix{Float64}(representation.coefficient_matrix),
        )
    elseif representation.contraction_kind == :identity
        return (
            states = _cartesian_product_states(representation.metadata.parent_axis_counts),
            coefficients = nothing,
        )
    end
    throw(
        ArgumentError(
            "unsupported Cartesian representation contraction kind :$(representation.contraction_kind)",
        ),
    )
end

function _cartesian_parent_cross_overlap(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    overlap_x = _basis_cross_overlap_1d(left.axis_representations.x, right.axis_representations.x)
    overlap_y = _basis_cross_overlap_1d(left.axis_representations.y, right.axis_representations.y)
    overlap_z = _basis_cross_overlap_1d(left.axis_representations.z, right.axis_representations.z)

    if left.contraction_kind == :identity && right.contraction_kind == :identity
        return Matrix{Float64}(kron(overlap_x, kron(overlap_y, overlap_z)))
    end

    left_parent = _cartesian_parent_state_basis(left)
    right_parent = _cartesian_parent_state_basis(right)
    raw_cross = _cartesian_product_cross_overlap(
        left_parent.states,
        right_parent.states,
        overlap_x,
        overlap_y,
        overlap_z,
    )

    if left_parent.coefficients !== nothing
        raw_cross = Matrix{Float64}(transpose(left_parent.coefficients) * raw_cross)
    end
    if right_parent.coefficients !== nothing
        raw_cross = Matrix{Float64}(raw_cross * right_parent.coefficients)
    end
    return raw_cross
end

"""
    cross_overlap(left::CartesianBasisRepresentation3D, right::CartesianBasisRepresentation3D)

Return the exact Cartesian basis cross-overlap matrix between `left` and
`right`, with rows ordered by the left basis and columns ordered by the right
basis.

This first pass supports representations whose defining parent/raw space is an
explicit Cartesian product basis. That covers:

- direct-product Cartesian representations
- nested fixed-block Cartesian representations
- exact cross overlaps between those families

The implementation is algebraic:

- if both sides expose an explicit Cartesian product parent, build the parent
  cross overlap from the stored 1D representation layers
- then contract by the stored 3D coefficient matrices

Hybrid/QW residual final bases are intentionally rejected for now. Their
current public representation does not yet carry enough exact raw
Cartesian-supplement cross-overlap metadata to make a fully honest exact path.
"""
function cross_overlap(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    _cartesian_supported_exact_parent_overlap(left, right) || throw(
        _cartesian_cross_overlap_error(left, right),
    )
    return _cartesian_parent_cross_overlap(left, right)
end

"""
    CartesianBasisTransferDiagnostics

Small public diagnostic bundle for exact Cartesian basis-to-basis transfer on
the supported pure Cartesian representation families.

The diagnostics record:

- the source and target basis dimensions
- whether the exact transfer used the same-parent or general pure-Cartesian
  path
- the residual of the Galerkin normal equations defining the transfer matrix
- symmetry defects of the source and target self-overlap matrices
- optional source/target metric traces and transfer residuals when specific
  orbital coefficients were transferred
"""
struct CartesianBasisTransferDiagnostics
    source_dimension::Int
    target_dimension::Int
    transfer_path::Symbol
    projector_residual_inf::Float64
    source_self_overlap_error::Float64
    target_self_overlap_error::Float64
    source_metric_trace::Union{Nothing,Float64}
    target_metric_trace::Union{Nothing,Float64}
    transferred_residual_inf::Union{Nothing,Float64}
end

function Base.show(io::IO, diagnostics::CartesianBasisTransferDiagnostics)
    print(
        io,
        "CartesianBasisTransferDiagnostics(source_dim=",
        diagnostics.source_dimension,
        ", target_dim=",
        diagnostics.target_dimension,
        ", path=:",
        diagnostics.transfer_path,
        ", projector_residual_inf=",
        diagnostics.projector_residual_inf,
        ")",
    )
end

"""
    CartesianBasisProjector3D

Exact metric-consistent Cartesian basis transfer matrix together with light
diagnostics describing how it was built.
"""
struct CartesianBasisProjector3D
    matrix::Matrix{Float64}
    diagnostics::CartesianBasisTransferDiagnostics
end

function Base.show(io::IO, projector::CartesianBasisProjector3D)
    print(
        io,
        "CartesianBasisProjector3D(size=",
        size(projector.matrix),
        ", path=:",
        projector.diagnostics.transfer_path,
        ")",
    )
end

"""
    CartesianOrbitalTransferResult

Returned by `transfer_orbitals(...)`. Bundles the transferred coefficients, the
exact basis projector used, and transfer diagnostics.
"""
struct CartesianOrbitalTransferResult{CT <: Union{Vector{Float64},Matrix{Float64}}}
    coefficients::CT
    projector::CartesianBasisProjector3D
    diagnostics::CartesianBasisTransferDiagnostics
end

function Base.show(io::IO, result::CartesianOrbitalTransferResult)
    print(
        io,
        "CartesianOrbitalTransferResult(size=",
        size(result.coefficients),
        ", path=:",
        result.diagnostics.transfer_path,
        ")",
    )
end

function _cartesian_transfer_path(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    return _cartesian_same_parent_raw_identity(source, target) ?
           :same_parent_metric_projection :
           :pure_cartesian_metric_projection
end

function _cartesian_transfer_diagnostics(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
    target_overlap::AbstractMatrix{<:Real},
    target_source_overlap::AbstractMatrix{<:Real},
    projector::AbstractMatrix{<:Real};
    source_coefficients::Union{Nothing,AbstractVecOrMat{<:Real}} = nothing,
    transferred_coefficients::Union{Nothing,AbstractVecOrMat{<:Real}} = nothing,
)
    source_overlap = cross_overlap(source, source)
    source_dimension = source.metadata.final_dimension
    target_dimension = target.metadata.final_dimension
    transfer_path = _cartesian_transfer_path(source, target)
    projector_residual_inf = norm(target_overlap * projector - target_source_overlap, Inf)
    source_self_overlap_error = norm(source_overlap - transpose(source_overlap), Inf)
    target_self_overlap_error = norm(target_overlap - transpose(target_overlap), Inf)

    if source_coefficients === nothing || transferred_coefficients === nothing
        return CartesianBasisTransferDiagnostics(
            source_dimension,
            target_dimension,
            transfer_path,
            projector_residual_inf,
            source_self_overlap_error,
            target_self_overlap_error,
            nothing,
            nothing,
            nothing,
        )
    end

    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_coefficients_matrix = _cartesian_coefficients_matrix(transferred_coefficients)
    source_metric =
        Matrix{Float64}(transpose(source_coefficients_matrix) * source_overlap * source_coefficients_matrix)
    target_metric =
        Matrix{Float64}(transpose(transferred_coefficients_matrix) * target_overlap * transferred_coefficients_matrix)
    transferred_residual_inf = norm(
        target_overlap * transferred_coefficients_matrix -
        target_source_overlap * source_coefficients_matrix,
        Inf,
    )

    return CartesianBasisTransferDiagnostics(
        source_dimension,
        target_dimension,
        transfer_path,
        projector_residual_inf,
        source_self_overlap_error,
        target_self_overlap_error,
        tr(source_metric),
        tr(target_metric),
        transferred_residual_inf,
    )
end

function _cartesian_coefficients_matrix(source_coefficients::AbstractVector{<:Real})
    return reshape(Float64[Float64(value) for value in source_coefficients], :, 1)
end

function _cartesian_coefficients_matrix(source_coefficients::AbstractMatrix{<:Real})
    return Matrix{Float64}(source_coefficients)
end

"""
    basis_projector(
        source::CartesianBasisRepresentation3D,
        target::CartesianBasisRepresentation3D,
    )

Return the exact metric-consistent basis-to-basis transfer operator from
`source` into `target` for the currently supported pure Cartesian
representations.

If `SAA = <A|A>`, `SBB = <B|B>`, and `SBA = <B|A>`, the projector matrix is the
Galerkin transfer

`T_{B <- A} = SBB \\ SBA`

so that target-space coefficients satisfy the normal equations

`SBB * C_B = SBA * C_A`.

This first pass supports the same pure Cartesian representation families as
`cross_overlap(...)`. Hybrid/QW residual final bases remain intentionally
unsupported here.
"""
function basis_projector(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    target_overlap = cross_overlap(target, target)
    target_source_overlap = cross_overlap(target, source)
    factorization = cholesky(Symmetric(Matrix{Float64}(target_overlap)))
    projector_matrix = Matrix{Float64}(factorization \ Matrix{Float64}(target_source_overlap))
    diagnostics = _cartesian_transfer_diagnostics(
        source,
        target,
        target_overlap,
        target_source_overlap,
        projector_matrix,
    )
    return CartesianBasisProjector3D(projector_matrix, diagnostics)
end

function _cartesian_transfer_coefficients(
    source_coefficients::AbstractVector{<:Real},
    projector::CartesianBasisProjector3D,
)
    source_dimension = size(projector.matrix, 2)
    length(source_coefficients) == source_dimension || throw(
        DimensionMismatch(
            "source coefficient vector length $(length(source_coefficients)) does not match source basis dimension $source_dimension",
        ),
    )
    return Vector{Float64}(projector.matrix * Vector{Float64}(source_coefficients))
end

function _cartesian_transfer_coefficients(
    source_coefficients::AbstractMatrix{<:Real},
    projector::CartesianBasisProjector3D,
)
    source_dimension = size(projector.matrix, 2)
    size(source_coefficients, 1) == source_dimension || throw(
        DimensionMismatch(
            "source coefficient matrix row count $(size(source_coefficients, 1)) does not match source basis dimension $source_dimension",
        ),
    )
    return Matrix{Float64}(projector.matrix * Matrix{Float64}(source_coefficients))
end

"""
    transfer_orbitals(
        source_coefficients::AbstractVecOrMat,
        source::CartesianBasisRepresentation3D,
        target::CartesianBasisRepresentation3D,
    )

Transfer orbital coefficients from `source` into `target` using the exact
metric-consistent projector returned by `basis_projector(...)`.

The returned coefficients satisfy the Galerkin normal equations

`SBB * C_B = SBA * C_A`

for the supported pure Cartesian representation families.
"""
function transfer_orbitals(
    source_coefficients::AbstractVecOrMat{<:Real},
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    projector = basis_projector(source, target)
    transferred_coefficients = _cartesian_transfer_coefficients(source_coefficients, projector)
    target_overlap = cross_overlap(target, target)
    target_source_overlap = cross_overlap(target, source)
    diagnostics = _cartesian_transfer_diagnostics(
        source,
        target,
        target_overlap,
        target_source_overlap,
        projector.matrix;
        source_coefficients = source_coefficients,
        transferred_coefficients = transferred_coefficients,
    )
    return CartesianOrbitalTransferResult(
        transferred_coefficients,
        projector,
        diagnostics,
    )
end

function _cartesian_basis_labels(prefix::AbstractString, count::Int)
    return [string(prefix, index) for index in 1:count]
end

function _cartesian_product_labels(
    orbitals::AbstractVector{<:CartesianProductOrbital3D},
)
    return [string("g(", orbital.ix, ",", orbital.iy, ",", orbital.iz, ")") for orbital in orbitals]
end

function _cartesian_labels(orbitals::AbstractVector{<:QiuWhiteHybridOrbital3D})
    return String[String(orbital.label) for orbital in orbitals]
end

function _cartesian_labels(orbitals::AbstractVector{<:_AtomicCartesianShellOrbital3D})
    return String[String(orbital.label) for orbital in orbitals]
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:CartesianProductOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.x
        centers_out[row, 2] = orbital.y
        centers_out[row, 3] = orbital.z
    end
    return centers_out
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:QiuWhiteHybridOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.x
        centers_out[row, 2] = orbital.y
        centers_out[row, 3] = orbital.z
    end
    return centers_out
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:_AtomicCartesianShellOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.center[1]
        centers_out[row, 2] = orbital.center[2]
        centers_out[row, 3] = orbital.center[3]
    end
    return centers_out
end

function _cartesian_axis_representations(basis::MappedUniformBasis)
    representation = basis_representation(basis)
    return (x = representation, y = representation, z = representation)
end

function _cartesian_axis_representations(basis::AbstractBondAlignedOrdinaryQWBasis3D)
    return (
        x = basis_representation(basis.basis_x),
        y = basis_representation(basis.basis_y),
        z = basis_representation(basis.basis_z),
    )
end

function _cartesian_axis_metadata(axis_representations::NamedTuple)
    return (
        x = axis_representations.x.metadata,
        y = axis_representations.y.metadata,
        z = axis_representations.z.metadata,
    )
end

function _cartesian_axis_signature(representation::BasisRepresentation1D)
    metadata = representation.metadata
    return (
        metadata.basis_kind,
        metadata.family_name,
        metadata.center_data,
        metadata.reference_center_data,
        metadata.integral_weight_data,
        metadata.basis_labels,
        metadata.coefficient_matrix,
    )
end

function _cartesian_axis_sharing(axis_representations::NamedTuple)
    x_signature = _cartesian_axis_signature(axis_representations.x)
    y_signature = _cartesian_axis_signature(axis_representations.y)
    z_signature = _cartesian_axis_signature(axis_representations.z)
    xy = isequal(x_signature, y_signature)
    xz = isequal(x_signature, z_signature)
    yz = isequal(y_signature, z_signature)
    xy && yz && return :shared_xyz
    xy && return :shared_xy
    xz && return :shared_xz
    yz && return :shared_yz
    return :separate_axes
end

function _cartesian_axis_counts(axis_representations::NamedTuple)
    return (
        length(axis_representations.x.metadata.center_data),
        length(axis_representations.y.metadata.center_data),
        length(axis_representations.z.metadata.center_data),
    )
end

function _cartesian_axis_centers(axis_representations::NamedTuple)
    return (
        axis_representations.x.metadata.center_data,
        axis_representations.y.metadata.center_data,
        axis_representations.z.metadata.center_data,
    )
end

function _cartesian_direct_product_representation(
    axis_representations::NamedTuple;
    route_metadata::NamedTuple = (;),
)
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    axis_counts = _cartesian_axis_counts(axis_representations)
    x_centers, y_centers, z_centers = _cartesian_axis_centers(axis_representations)
    orbitals = _mapped_cartesian_orbitals(x_centers, y_centers, z_centers)
    labels = _cartesian_product_labels(orbitals)
    centers_matrix = _cartesian_center_matrix(orbitals)
    metadata = CartesianBasisMetadata3D(
        :direct_product,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_product_basis,
        axis_counts,
        length(orbitals),
        length(orbitals),
        nothing,
        labels,
        centers_matrix,
        route_metadata,
    )
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :identity,
        nothing,
        copy(labels),
        copy(centers_matrix),
        nothing,
        nothing,
        (;),
    )
end

function _cartesian_direct_product_representation(
    basis::MappedUniformBasis;
    route_metadata::NamedTuple = (;),
)
    return _cartesian_direct_product_representation(
        _cartesian_axis_representations(basis);
        route_metadata = route_metadata,
    )
end

function _cartesian_direct_product_representation(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    route_metadata::NamedTuple = (;),
)
    return _cartesian_direct_product_representation(
        _cartesian_axis_representations(basis);
        route_metadata = route_metadata,
    )
end

function _cartesian_basis_route_metadata(basis::BondAlignedDiatomicQWBasis3D)
    return (
        basis_family = :bond_aligned_diatomic,
        bond_axis = basis.bond_axis,
        nuclei = copy(basis.nuclei),
        target_core_spacing = basis.target_core_spacing,
    )
end

function _cartesian_basis_route_metadata(basis::BondAlignedHomonuclearChainQWBasis3D)
    return (
        basis_family = :bond_aligned_homonuclear_chain,
        chain_axis = basis.chain_axis,
        chain_coordinates = Float64[Float64(value) for value in basis.chain_coordinates],
        nuclei = copy(basis.nuclei),
        target_core_spacing = basis.target_core_spacing,
    )
end

function _cartesian_basis_route_metadata(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    return (
        basis_family = :axis_aligned_homonuclear_square_lattice,
        lattice_size = basis.lattice_size,
        x_coordinates = Float64[Float64(value) for value in basis.x_coordinates],
        y_coordinates = Float64[Float64(value) for value in basis.y_coordinates],
        nuclei = copy(basis.nuclei),
        target_core_spacing = basis.target_core_spacing,
    )
end

function basis_representation(basis::BondAlignedDiatomicQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

function basis_representation(basis::BondAlignedHomonuclearChainQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

function basis_representation(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

basis_metadata(basis::AbstractBondAlignedOrdinaryQWBasis3D) = basis_representation(basis).metadata

function _cartesian_parent_axis_representations(parent_basis)
    parent_basis isa MappedUniformBasis && return _cartesian_axis_representations(parent_basis)
    parent_basis isa AbstractBondAlignedOrdinaryQWBasis3D &&
        return _cartesian_axis_representations(parent_basis)
    throw(
        ArgumentError(
            "Cartesian basis representation does not yet support parent basis $(typeof(parent_basis))",
        ),
    )
end

function _cartesian_full_parent_box(axis_counts::NTuple{3,Int})
    return (1:axis_counts[1], 1:axis_counts[2], 1:axis_counts[3])
end

function _cartesian_shell_kind(::Any)
    return :unknown_shell
end

function _cartesian_shell_kind(::_CartesianNestedShell3D)
    return :shell
end

function _cartesian_shell_kind(::_CartesianNestedCompleteShell3D)
    return :complete_shell
end

function _cartesian_shell_kind(::_CartesianNestedShellPlusCore3D)
    return :shell_plus_core
end

function _cartesian_shell_kind(::_CartesianNestedShellSequence3D)
    return :shell_sequence
end

function _cartesian_shell_working_box(shell)
    return hasproperty(shell, :working_box) ? shell.working_box : nothing
end

function _cartesian_working_box_profile(
    working_box::Union{Nothing,NTuple{3,UnitRange{Int}}},
    axis_counts::NTuple{3,Int},
)
    working_box === nothing && return nothing
    return working_box == _cartesian_full_parent_box(axis_counts) ? :full_parent : :explicit_inner_box
end

function _cartesian_maybe_complete_shell_nside(shell)
    shell isa _CartesianNestedShellSequence3D || return nothing
    working_box = shell.working_box
    length(working_box[1]) == length(working_box[2]) == length(working_box[3]) || return nothing
    core_side = round(Int, cbrt(length(shell.core_indices)))
    core_side^3 == length(shell.core_indices) || return nothing
    expected_increment = _one_center_atomic_shell_increment(core_side)
    all(length(range) == expected_increment for range in shell.layer_column_ranges) || return nothing
    return core_side
end

function _cartesian_support_states(
    fixed_block::_NestedFixedBlock3D,
    axis_counts::NTuple{3,Int},
)
    if hasproperty(fixed_block.shell, :support_states)
        return NTuple{3,Int}[state for state in fixed_block.shell.support_states]
    end
    return NTuple{3,Int}[
        _cartesian_unflat_index(index, axis_counts) for index in fixed_block.support_indices
    ]
end

function _cartesian_fixed_block_route_metadata(
    fixed_block::_NestedFixedBlock3D,
    axis_counts::NTuple{3,Int},
)
    working_box = _cartesian_shell_working_box(fixed_block.shell)
    return (
        shell_kind = _cartesian_shell_kind(fixed_block.shell),
        working_box_profile = _cartesian_working_box_profile(working_box, axis_counts),
        nside = _cartesian_maybe_complete_shell_nside(fixed_block.shell),
        support_count = length(fixed_block.support_indices),
        term_storage = fixed_block.term_storage,
    )
end

function basis_representation(fixed_block::_NestedFixedBlock3D)
    axis_representations = _cartesian_parent_axis_representations(fixed_block.parent_basis)
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    axis_counts = _cartesian_axis_counts(axis_representations)
    x_centers, y_centers, z_centers = _cartesian_axis_centers(axis_representations)
    parent_orbitals = _mapped_cartesian_orbitals(x_centers, y_centers, z_centers)
    parent_labels = _cartesian_product_labels(parent_orbitals)
    parent_centers = _cartesian_center_matrix(parent_orbitals)
    support_states = _cartesian_support_states(fixed_block, axis_counts)
    final_labels = _cartesian_basis_labels("nf", size(fixed_block.coefficient_matrix, 2))
    metadata = CartesianBasisMetadata3D(
        :nested_fixed_block,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_product_basis,
        axis_counts,
        length(parent_orbitals),
        size(fixed_block.coefficient_matrix, 2),
        _cartesian_shell_working_box(fixed_block.shell),
        final_labels,
        Matrix{Float64}(fixed_block.fixed_centers),
        _cartesian_fixed_block_route_metadata(fixed_block, axis_counts),
    )
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :dense,
        Matrix{Float64}(fixed_block.coefficient_matrix),
        parent_labels,
        parent_centers,
        Vector{Int}(fixed_block.support_indices),
        support_states,
        (;),
    )
end

basis_metadata(fixed_block::_NestedFixedBlock3D) = basis_representation(fixed_block).metadata

function _cartesian_empty_centers()
    return zeros(Float64, 0, 3)
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
        labels = String[],
        centers = _cartesian_empty_centers(),
        angular_powers = NTuple{3,Int}[],
        exponents = Vector{Vector{Float64}}(),
        coefficients = Vector{Vector{Float64}}(),
    )
end

function _cartesian_supplement_metadata(
    data::Union{
        LegacyAtomicGaussianSupplement,
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    supplement3d =
        data isa LegacyAtomicGaussianSupplement ? _atomic_cartesian_shell_supplement_3d(data) :
        _bond_aligned_diatomic_cartesian_shell_supplement_3d(data)
    orbitals = supplement3d.orbitals
    return (
        labels = _cartesian_labels(orbitals),
        centers = _cartesian_center_matrix(orbitals),
        angular_powers = NTuple{3,Int}[(orbital.lx, orbital.ly, orbital.lz) for orbital in orbitals],
        exponents = [Float64[orbital.exponents...] for orbital in orbitals],
        coefficients = [Float64[orbital.coefficients...] for orbital in orbitals],
    )
end

function _qwrg_cartesian_parent_representation(
    operators::QiuWhiteResidualGaussianOperators,
)
    basis = operators.basis
    basis isa _NestedFixedBlock3D && return basis_representation(basis)
    basis isa AbstractBondAlignedOrdinaryQWBasis3D && return basis_representation(basis)
    basis isa MappedUniformBasis && return _cartesian_direct_product_representation(basis)
    throw(
        ArgumentError(
            "Cartesian basis representation does not yet support QW parent basis $(typeof(basis))",
        ),
    )
end

function basis_representation(operators::QiuWhiteResidualGaussianOperators)
    cartesian_parent = _qwrg_cartesian_parent_representation(operators)
    axis_representations = cartesian_parent.axis_representations
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    size(cartesian_parent.metadata.basis_centers, 1) == operators.gausslet_count || throw(
        ArgumentError(
            "QW Cartesian basis representation requires the parent Cartesian representation to match gausslet_count",
        ),
    )
    supplement_metadata = _cartesian_supplement_metadata(operators.gaussian_data)
    parent_labels = vcat(
        cartesian_parent.metadata.basis_labels,
        supplement_metadata.labels,
    )
    parent_centers = vcat(
        cartesian_parent.metadata.basis_centers,
        supplement_metadata.centers,
    )
    size(operators.raw_to_final, 1) == length(parent_labels) || throw(
        ArgumentError(
            "QW Cartesian basis representation requires raw_to_final rows to match the explicit raw parent basis dimension",
        ),
    )
    final_labels = _cartesian_labels(operators.orbital_data)
    final_centers = _cartesian_center_matrix(operators.orbital_data)
    metadata = CartesianBasisMetadata3D(
        :hybrid_residual,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_plus_supplement_raw,
        cartesian_parent.metadata.parent_axis_counts,
        size(operators.raw_to_final, 1),
        size(operators.raw_to_final, 2),
        cartesian_parent.metadata.working_box,
        final_labels,
        final_centers,
        (
            gausslet_count = operators.gausslet_count,
            residual_count = operators.residual_count,
            supplement_kind = _cartesian_supplement_kind(operators.gaussian_data),
            supplement_lmax = _cartesian_supplement_lmax(operators.gaussian_data),
        ),
    )
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :dense,
        Matrix{Float64}(operators.raw_to_final),
        parent_labels,
        parent_centers,
        nothing,
        nothing,
        (
            cartesian_parent_representation = cartesian_parent,
            supplement_orbitals = supplement_metadata,
        ),
    )
end

basis_metadata(operators::QiuWhiteResidualGaussianOperators) = basis_representation(operators).metadata
