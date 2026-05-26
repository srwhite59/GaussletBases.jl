using Printf
using LinearAlgebra
using SparseArrays

using GaussletBases

const CP = GaussletBases.CartesianParentGaussletBases
const CCP = GaussletBases.CartesianContractedParents
const CCPM = GaussletBases.CartesianContractedParentMetrics

function _synthetic_sparse_coefficients(nparent::Int, ncols::Int, nnz_per_col::Int)
    rows = Int[]
    cols = Int[]
    values = Float64[]
    sizehint!(rows, ncols * nnz_per_col)
    sizehint!(cols, ncols * nnz_per_col)
    sizehint!(values, ncols * nnz_per_col)
    for col in 1:ncols, entry in 1:nnz_per_col
        row = mod(97 * col + 389 * entry + 17 * col * entry, nparent) + 1
        value = (-1.0)^(entry + col) / sqrt(float(nnz_per_col))
        push!(rows, row)
        push!(cols, col)
        push!(values, value)
    end
    return sparse(rows, cols, values, nparent, ncols)
end

axis = build_basis(
    MappedUniformBasisSpec(
        :G10;
        count = 15,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ),
)
parent = CP.cartesian_parent_gausslet_basis(axis)
parent_dim = CP.parent_dimension(parent)
contracted_dim = 96
nnz_per_col = 10
coefficients = _synthetic_sparse_coefficients(parent_dim, contracted_dim, nnz_per_col)
contracted = CCP.CartesianContractedParent3D(
    parent,
    coefficients;
    metadata = (source = :synthetic_metric_packet_benchmark,),
)

# Compile first, then report the steady small-fixture timing.
warm_packet = CCPM.cartesian_contracted_parent_metric_packet(contracted)
timed = @timed CCPM.cartesian_contracted_parent_metric_packet(contracted)
packet = timed.value

println("Cartesian contracted parent metric packet benchmark")
@printf("parent_side = %d\n", 15)
@printf("parent_dimension = %d\n", parent_dim)
@printf("contracted_dimension = %d\n", contracted_dim)
@printf("coefficient_nnz = %d\n", nnz(coefficients))
@printf("max_column_nnz = %d\n", packet.diagnostics.max_column_nnz)
@printf("construction_path = %s\n", String(packet.diagnostics.construction_path))
@printf("dense_parent_matrix_used = %s\n", string(packet.diagnostics.dense_parent_matrix_used))
@printf("overlap_symmetry_error = %.6e\n", packet.diagnostics.overlap_symmetry_error)
@printf("overlap_identity_error = %.6e\n", packet.diagnostics.overlap_identity_error)
@printf("warm_time_s = %.6f\n", timed.time)
@printf("warm_alloc_mib = %.3f\n", timed.bytes / 1024^2)
@printf("warmup_overlap_trace = %.12e\n", tr(warm_packet.overlap))
