# Integration/slow test. Do not include in default nested runner.

@testset "Nested coefficient maps support sparse storage" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    direct_fixed_block = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    sparse_fixed_block = _with_sparse_nested_coefficients(direct_fixed_block)

    direct_representation = basis_representation(direct_fixed_block)
    sparse_representation = basis_representation(sparse_fixed_block)
    support_coefficients = GaussletBases._nested_support_coefficient_slice(
        direct_fixed_block.shell.coefficient_matrix,
        direct_fixed_block.shell.support_indices,
    )
    support_workspace, contraction_scratch = GaussletBases._nested_support_reference_workspaces(
        support_coefficients,
        length(direct_fixed_block.shell.support_indices),
        size(direct_fixed_block.shell.coefficient_matrix, 2),
    )

    @test direct_fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test direct_representation.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sparse_fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sparse_representation.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test support_coefficients isa SparseMatrixCSC{Float64,Int}
    @test size(support_workspace) == (0, 0)
    @test size(contraction_scratch) == (0, 0)
    @test size(support_coefficients) == (
        length(direct_fixed_block.shell.support_indices),
        size(direct_fixed_block.shell.coefficient_matrix, 2),
    )
    @test Matrix(sparse_representation.coefficient_matrix) ≈
        Matrix(direct_representation.coefficient_matrix) atol = 1.0e-12 rtol = 1.0e-12
    @test cross_overlap(sparse_representation, sparse_representation) ≈
        cross_overlap(direct_representation, direct_representation) atol = 1.0e-10 rtol = 1.0e-10

    mktemp() do sparse_path, sparse_io
        close(sparse_io)
        sparse_matrix = sparse(direct_fixed_block.coefficient_matrix)
        jldopen(sparse_path, "w") do file
            file["matrix"] = sparse_matrix
        end
        restored = jldopen(sparse_path, "r") do file
            file["matrix"]
        end
        @test restored isa SparseMatrixCSC{Float64,Int}
        @test restored == sparse_matrix
    end
end

function _atomic_hybrid_cartesian_representation_fixture()
    return _cached_fixture(:atomic_hybrid_cartesian_representation_fixture, () -> begin
        basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = 13,
                mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
                reference_spacing = 1.0,
            ),
        )
        expansion = coulomb_gaussian_expansion(doacc = false)
        fixed_full = one_center_atomic_full_parent_fixed_block(
            basis;
            expansion,
            nside = 5,
        )
        fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
            basis;
            expansion,
            working_box = 2:12,
            nside = 5,
        )
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        full_ops = ordinary_cartesian_qiu_white_operators(
            fixed_full,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        legacy_ops = ordinary_cartesian_qiu_white_operators(
            fixed_legacy,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        (
            basis = basis,
            expansion = expansion,
            fixed_full = fixed_full,
            fixed_legacy = fixed_legacy,
            fixed_full_rep = basis_representation(fixed_full),
            fixed_legacy_rep = basis_representation(fixed_legacy),
            supplement = supplement,
            full_ops = full_ops,
            legacy_ops = legacy_ops,
            full_rep = basis_representation(full_ops),
            legacy_rep = basis_representation(legacy_ops),
        )
    end)
end

function _metric_normalize_orbital(
    coefficients::AbstractVector,
    overlap::AbstractMatrix{<:Real},
)
    orbital = Float64[Float64(real(value)) for value in coefficients]
    norm2 = Float64(real(dot(orbital, overlap * orbital)))
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero target-metric norm"))
    return orbital ./ sqrt(norm2)
end

function _metric_orbital_overlap(
    left::AbstractVector,
    right::AbstractVector,
    overlap::AbstractMatrix{<:Real},
)
    normalized_left = _metric_normalize_orbital(left, overlap)
    normalized_right = _metric_normalize_orbital(right, overlap)
    return Float64(real(dot(normalized_left, overlap * normalized_right)))
end

function _ordinary_cartesian_hybrid_orbital_observables(
    operators::OrdinaryCartesianOperators3D,
    orbital::AbstractVector;
    overlap_tol::Real = 1.0e-7,
)
    overlap = Matrix{Float64}(operators.overlap)
    normalized = _metric_normalize_orbital(orbital, overlap)
    one_body = Float64(real(dot(normalized, operators.one_body_hamiltonian * normalized)))
    vee = GaussletBases.ordinary_cartesian_vee_expectation(
        operators,
        normalized;
        overlap_tol = overlap_tol,
    )
    return (
        orbital = normalized,
        metric_norm_error = abs(Float64(real(dot(normalized, overlap * normalized))) - 1.0),
        one_body = one_body,
        vee = vee,
        total = 2.0 * one_body + vee,
    )
end

function _atomic_direct_product_he_extent_change_contract_fixture(;
    source_count::Int = 3,
    target_count::Int = 5,
)
    key = Symbol(:atomic_direct_product_he_extent_change_contract, source_count, target_count)
    return _cached_fixture(key, () -> begin
        mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
        source_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = source_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )
        target_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = target_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )

        source_rep = basis_representation(source_basis)
        target_rep = basis_representation(target_basis)
        offset = (target_count - source_count) ÷ 2
        shared_slice = (offset + 1):(offset + source_count)

        return (
            source_count = source_count,
            target_count = target_count,
            shared_slice = shared_slice,
            source_rep = source_rep,
            target_rep = target_rep,
            centers_subset =
                source_rep.metadata.center_data == target_rep.metadata.center_data[shared_slice],
            weights_subset =
                source_rep.metadata.integral_weight_data ==
                target_rep.metadata.integral_weight_data[shared_slice],
            coefficient_core_match =
                source_rep.coefficient_matrix ==
                target_rep.coefficient_matrix[shared_slice, shared_slice],
        )
    end)
end

function _atomic_hybrid_he_same_parent_stress_fixture(;
    parent_count::Int = 7,
    source_working_box::UnitRange{Int} = 2:6,
    supplement_lmax::Int = 1,
)
    key = Symbol(
        :atomic_hybrid_he_same_parent_stress_fixture,
        parent_count,
        first(source_working_box),
        last(source_working_box),
        supplement_lmax,
    )
    return _cached_fixture(key, () -> begin
        expansion = coulomb_gaussian_expansion(doacc = false)
        mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = supplement_lmax)

        parent_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = parent_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )

        source_fixed = one_center_atomic_legacy_profile_fixed_block(
            parent_basis;
            expansion,
            working_box = source_working_box,
            nside = 5,
        )
        target_fixed = one_center_atomic_full_parent_fixed_block(
            parent_basis;
            expansion,
            nside = 5,
        )

        source_ops = ordinary_cartesian_qiu_white_operators(
            source_fixed,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        target_ops = ordinary_cartesian_qiu_white_operators(
            target_fixed,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )

        source_rep = basis_representation(source_ops)
        target_rep = basis_representation(target_ops)
        source_check = GaussletBases.ordinary_cartesian_1s2_check(
            source_ops;
            overlap_tol = 1.0e-7,
        )
        target_check = GaussletBases.ordinary_cartesian_1s2_check(
            target_ops;
            overlap_tol = 1.0e-7,
        )
        source_observables = _ordinary_cartesian_hybrid_orbital_observables(
            source_ops,
            source_check.orbital;
            overlap_tol = 1.0e-7,
        )
        target_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            target_check.orbital;
            overlap_tol = 1.0e-7,
        )

        transfer = transfer_orbitals(source_observables.orbital, source_rep, target_rep)
        transferred_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            transfer.coefficients;
            overlap_tol = 1.0e-7,
        )
        target_overlap = Matrix{Float64}(target_ops.overlap)
        overlap_with_target = _metric_orbital_overlap(
            transferred_observables.orbital,
            target_observables.orbital,
            target_overlap,
        )
        sign = overlap_with_target < 0.0 ? -1.0 : 1.0
        aligned_transferred_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            sign .* transferred_observables.orbital;
            overlap_tol = 1.0e-7,
        )
        aligned_overlap_to_target = abs(
            _metric_orbital_overlap(
                aligned_transferred_observables.orbital,
                target_observables.orbital,
                target_overlap,
            ),
        )

        return (
            parent_basis = parent_basis,
            source_fixed = source_fixed,
            target_fixed = target_fixed,
            supplement = supplement,
            source_working_box = source_working_box,
            target_working_box = target_fixed.shell.working_box,
            source_ops = source_ops,
            target_ops = target_ops,
            source_rep = source_rep,
            target_rep = target_rep,
            source_check = source_check,
            target_check = target_check,
            source_observables = source_observables,
            target_observables = target_observables,
            transfer = transfer,
            transferred_observables = transferred_observables,
            aligned_transferred_observables = aligned_transferred_observables,
            aligned_overlap_to_target = aligned_overlap_to_target,
        )
    end)
end
