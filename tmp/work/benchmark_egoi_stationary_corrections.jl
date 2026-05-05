using LinearAlgebra
using Printf

using GaussletBases

function _symmetric_reference(n::Int)
    matrix = [sin(0.17 * row + 0.31 * column) for row in 1:n, column in 1:n]
    return 0.5 .* (matrix .+ transpose(matrix))
end

function main()
    nbasis = 40
    ntarget = 3
    qraw = [cos(0.11 * row * column) + sin(0.07 * row + column) for row in 1:nbasis, column in 1:ntarget]
    qtarget = Matrix(qr(qraw).Q)[:, 1:ntarget]
    h = _symmetric_reference(nbasis)
    v = _symmetric_reference(nbasis)
    for index in 1:nbasis
        v[index, index] += 2.0
    end
    known_delta = 0.01 .* _symmetric_reference(nbasis)
    product = egoi_target_product_matrix(qtarget)
    exact_target = transpose(product) * (v + known_delta) * product
    occupations = [2.0, 2.0, 0.0]

    GC.gc()
    timed = @timed egoi_stationary_hamiltonian_correction(
        h,
        v,
        qtarget,
        exact_target,
        occupations,
    )
    result = timed.value

    println("case=synthetic_orthonormal_working_basis")
    println("timing_kind=fresh_process")
    println("nbasis=$(nbasis)")
    println("ntarget=$(ntarget)")
    println("target_pair_count=$(ntarget^2)")
    @printf("time_s=%.6f alloc_mb=%.3f gc_s=%.6f\n", timed.time, timed.bytes / 1024^2, timed.gctime)
    @printf("egoi_residual_max_after=%.6e\n", result.diagnostics.egoi.target_residual_max_after)
    @printf(
        "stationary_residual_max_after=%.6e\n",
        result.diagnostics.stationary.occupied_virtual_residual_max_after,
    )
    @printf("delta_v_fro=%.6e\n", norm(result.interaction_delta))
    @printf("delta_h_fro=%.6e\n", norm(result.one_body_delta))
end

main()
