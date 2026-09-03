"""
    CartesianBasisTransferDiagnostics

Small public diagnostic bundle for exact Cartesian basis-to-basis transfer on
the supported Cartesian representation families.

The diagnostics record:

- the source and target basis dimensions
- whether the exact transfer used the same-parent pure-Cartesian,
  general pure-Cartesian, or hybrid mixed-raw path
- the defect of the materialized transfer operator relative to the exact final
  basis cross overlap used by the working transfer contract
- optional self-overlap diagnostics when that data was explicitly constructed
- optional source/target Euclidean norm traces and transfer residuals when
  specific orbital coefficients were transferred
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

Returned by `transfer_orbitals(...)`. Bundles the transferred coefficients,
the exact basis projector when one was materialized, and transfer diagnostics.

`projector === nothing` means the caller selected a supported no-projector
transfer path, so the exact `C_B = S_BA * C_A` transfer was applied directly to
the supplied coefficient block without storing the dense `S_BA` matrix.
"""
struct CartesianOrbitalTransferResult{CT <: Union{Vector{Float64},Matrix{Float64}}}
    coefficients::CT
    projector::Union{Nothing,CartesianBasisProjector3D}
    diagnostics::CartesianBasisTransferDiagnostics
end

function Base.show(io::IO, result::CartesianOrbitalTransferResult)
    print(
        io,
        "CartesianOrbitalTransferResult(size=",
        size(result.coefficients),
        ", path=:",
        result.diagnostics.transfer_path,
        ", projector_size=",
        result.projector === nothing ? "none" : string(size(result.projector.matrix)),
        ")",
    )
end

function _cartesian_transfer_path(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    if source.metadata.parent_kind == :cartesian_plus_supplement_raw ||
       target.metadata.parent_kind == :cartesian_plus_supplement_raw
        return :hybrid_mixed_raw_cross_overlap_transfer
    end
    return _cartesian_same_parent_raw_identity(source, target) ?
           :same_parent_cross_overlap_transfer :
           :pure_cartesian_cross_overlap_transfer
end

function _cartesian_columns_are_nearly_orthonormal(
    coefficients::AbstractMatrix{<:Real};
    atol::Float64 = 1.0e-8,
)
    size(coefficients, 2) > 1 || return false
    gram = Matrix{Float64}(transpose(coefficients) * coefficients)
    identity = Matrix{Float64}(LinearAlgebra.I, size(gram, 1), size(gram, 2))
    return norm(gram - identity, Inf) <= atol
end

function _cartesian_euclidean_orthonormalize_columns(
    coefficients::AbstractMatrix{<:Real},
)
    ncols = size(coefficients, 2)
    factorization = qr(Matrix{Float64}(coefficients))
    orthonormal = Matrix{Float64}(factorization.Q[:, 1:ncols])
    rblock = Matrix{Float64}(factorization.R[1:ncols, 1:ncols])
    for column in 1:ncols
        if rblock[column, column] < 0.0
            orthonormal[:, column] .*= -1.0
        end
    end
    return orthonormal
end

function _cartesian_finalize_transferred_coefficients(
    source_coefficients::AbstractVector{<:Real},
    transferred_coefficients_matrix::AbstractMatrix{<:Real},
)
    return (
        vec(Matrix{Float64}(transferred_coefficients_matrix)),
        false,
    )
end

function _cartesian_finalize_transferred_coefficients(
    source_coefficients::AbstractMatrix{<:Real},
    transferred_coefficients_matrix::AbstractMatrix{<:Real},
)
    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_matrix = Matrix{Float64}(transferred_coefficients_matrix)
    if _cartesian_columns_are_nearly_orthonormal(source_coefficients_matrix) &&
       !_cartesian_columns_are_nearly_orthonormal(transferred_matrix)
        return (
            _cartesian_euclidean_orthonormalize_columns(transferred_matrix),
            true,
        )
    end
    return (transferred_matrix, false)
end

function _cartesian_transfer_diagnostics(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
    target_source_overlap::AbstractMatrix{<:Real},
    projector::AbstractMatrix{<:Real};
    source_coefficients::Union{Nothing,AbstractVecOrMat{<:Real}} = nothing,
    transferred_coefficients::Union{Nothing,AbstractVecOrMat{<:Real}} = nothing,
    orthonormalized::Bool = false,
)
    source_dimension = source.metadata.final_dimension
    target_dimension = target.metadata.final_dimension
    transfer_path = _cartesian_transfer_path(source, target)
    projector_residual_inf = norm(projector - target_source_overlap, Inf)

    if source_coefficients === nothing || transferred_coefficients === nothing
        return CartesianBasisTransferDiagnostics(
            source_dimension,
            target_dimension,
            transfer_path,
            projector_residual_inf,
            NaN,
            NaN,
            nothing,
            nothing,
            nothing,
        )
    end

    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_coefficients_matrix = _cartesian_coefficients_matrix(transferred_coefficients)
    source_metric =
        Matrix{Float64}(transpose(source_coefficients_matrix) * source_coefficients_matrix)
    target_metric =
        Matrix{Float64}(transpose(transferred_coefficients_matrix) * transferred_coefficients_matrix)
    transferred_residual_inf =
        orthonormalized ? NaN : norm(
            transferred_coefficients_matrix - projector * source_coefficients_matrix,
            Inf,
        )

    return CartesianBasisTransferDiagnostics(
        source_dimension,
        target_dimension,
        transfer_path,
        projector_residual_inf,
        NaN,
        NaN,
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

function _cartesian_transfer_diagnostics_without_projector(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
    source_coefficients::AbstractVecOrMat{<:Real},
    transferred_coefficients::AbstractVecOrMat{<:Real},
)
    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_coefficients_matrix = _cartesian_coefficients_matrix(transferred_coefficients)
    return CartesianBasisTransferDiagnostics(
        source.metadata.final_dimension,
        target.metadata.final_dimension,
        _cartesian_transfer_path(source, target),
        NaN,
        NaN,
        NaN,
        tr(transpose(source_coefficients_matrix) * source_coefficients_matrix),
        tr(transpose(transferred_coefficients_matrix) * transferred_coefficients_matrix),
        NaN,
    )
end

function _cartesian_basis_projector_with_stage_timings(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    overlap_result =
        _cartesian_uses_exact_mixed_raw_cross_overlap(target, source) ?
        _cartesian_mixed_raw_cross_overlap_with_stage_timings(target, source) :
        (matrix = Matrix{Float64}(cross_overlap(target, source)), stage_timings = nothing)
    diagnostics = _cartesian_transfer_diagnostics(
        source,
        target,
        overlap_result.matrix,
        overlap_result.matrix,
    )
    return (
        projector = CartesianBasisProjector3D(overlap_result.matrix, diagnostics),
        stage_timings = overlap_result.stage_timings,
    )
end

"""
    basis_projector(
        source::CartesianBasisRepresentation3D,
        target::CartesianBasisRepresentation3D,
    )

Return the exact metric-consistent basis-to-basis transfer operator from
`source` into `target` for the currently supported Cartesian representations.

For the final working-basis contract used in this repo, the transfer matrix is
the exact final-basis cross overlap

`T_{B <- A} = S_{BA} = <B|A>`.

This first pass supports the same representation families as `cross_overlap(...)`,
including hybrid residual-Gaussian final bases whose public representation
exposes the exact mixed raw-space identity.

Final-basis self-overlaps are construction/debug diagnostics only. Normal
transfer on these orthonormal working bases uses only the exact cross overlap;
raw/nonorthogonal generalized-overlap algebra remains outside this helper.
"""
function basis_projector(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    return _cartesian_basis_projector_with_stage_timings(source, target).projector
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

function transfer_orbitals(
    source_coefficients::AbstractVecOrMat{<:Real},
    projector::CartesianBasisProjector3D,
)
    transferred_coefficients_matrix = _cartesian_transfer_coefficients(source_coefficients, projector)
    transferred_coefficients, orthonormalized =
        _cartesian_finalize_transferred_coefficients(
            source_coefficients,
            _cartesian_coefficients_matrix(transferred_coefficients_matrix),
        )
    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_coefficients_matrix = _cartesian_coefficients_matrix(transferred_coefficients)
    diagnostics = CartesianBasisTransferDiagnostics(
        projector.diagnostics.source_dimension,
        projector.diagnostics.target_dimension,
        projector.diagnostics.transfer_path,
        projector.diagnostics.projector_residual_inf,
        projector.diagnostics.source_self_overlap_error,
        projector.diagnostics.target_self_overlap_error,
        tr(transpose(source_coefficients_matrix) * source_coefficients_matrix),
        tr(transpose(transferred_coefficients_matrix) * transferred_coefficients_matrix),
        orthonormalized ? NaN : norm(
            transferred_coefficients_matrix -
            projector.matrix * source_coefficients_matrix,
            Inf,
        ),
    )
    return CartesianOrbitalTransferResult(
        transferred_coefficients,
        projector,
        diagnostics,
    )
end

"""
    transfer_orbitals(
        source_coefficients::AbstractVecOrMat,
        source::CartesianBasisRepresentation3D,
        target::CartesianBasisRepresentation3D,
        ;
        materialize_projector::Bool = true,
    )

For the final orthonormal working-basis contract used in this repo, the transfer
formula is

`C_B = S_{BA} * C_A`.

With the default `materialize_projector = true`, this preserves the historical
explicit-projector route: build the exact dense projector `S_BA` with
`basis_projector(...)`, then apply it to `source_coefficients`.

With `materialize_projector = false`, supported mixed-raw hybrid transfers apply
the same exact formula directly to the supplied coefficient block and may return
`CartesianOrbitalTransferResult.projector === nothing`. Unsupported source/target
pairs throw rather than silently falling back to projector materialization.

Final-basis self-overlaps are diagnostic only and are not part of this normal
transfer path. When the source columns already represent an orthonormal occupied
block, the transferred target block may be Euclidean-orthonormalized as one
final cleanup step.
"""
function transfer_orbitals(
    source_coefficients::AbstractVecOrMat{<:Real},
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D;
    materialize_projector::Bool = true,
)
    if !materialize_projector
        _cartesian_uses_exact_mixed_raw_cross_overlap(target, source) || throw(
            ArgumentError(
                "materialize_projector=false currently requires an exact mixed-raw Cartesian transfer path",
            ),
        )
        source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
        transfer_result = _cartesian_mixed_raw_cross_apply_with_stage_timings(
            target,
            source,
            source_coefficients_matrix,
        )
        transferred_coefficients, _ =
            _cartesian_finalize_transferred_coefficients(
                source_coefficients,
                transfer_result.coefficients,
            )
        return CartesianOrbitalTransferResult(
            transferred_coefficients,
            nothing,
            _cartesian_transfer_diagnostics_without_projector(
                source,
                target,
                source_coefficients,
                transferred_coefficients,
            ),
        )
    end
    projector = basis_projector(source, target)
    return transfer_orbitals(source_coefficients, projector)
end
