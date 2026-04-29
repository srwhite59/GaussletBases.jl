using LinearAlgebra
using Statistics
using GaussletBases

BLAS.set_num_threads(1)

function _parse_args(args)
    values = Dict{String,String}()
    for arg in args
        key_value = split(arg, '='; limit = 2)
        length(key_value) == 2 || throw(
            ArgumentError("expected key=value argument, got $(repr(arg))"),
        )
        values[key_value[1]] = key_value[2]
    end
    return values
end

function _rss_bytes()
    isdefined(Sys, :maxrss) || return nothing
    try
        value = Sys.maxrss()
        value isa Integer || return nothing
        return Int(value)
    catch
        return nothing
    end
end

function _timed_stage(thunk)
    GC.gc()
    rss_before = _rss_bytes()
    result = nothing
    elapsed_seconds = @elapsed begin
        result = thunk()
    end
    rss_after = _rss_bytes()
    return (
        result = result,
        seconds = elapsed_seconds,
        rss_before_bytes = rss_before,
        rss_after_bytes = rss_after,
    )
end

function _coefficient_max_abs_residual(
    left::AbstractVecOrMat{<:Real},
    right::AbstractVecOrMat{<:Real},
)
    return maximum(abs.(Matrix{Float64}(left) .- Matrix{Float64}(right)))
end

function _case_paths()
    dir = joinpath(pwd(), "tmp", "hybrid_ns5_benchmark")
    return (
        dir = dir,
        source = joinpath(dir, "he_full_ns5.jld2"),
        target = joinpath(dir, "he_legacy_ns5.jld2"),
    )
end

function _he_ns5_fixture()
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
        expansion = expansion,
        nside = 5,
    )
    fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion = expansion,
        working_box = 2:12,
        nside = 5,
    )
    supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
    full_ops = ordinary_cartesian_qiu_white_operators(
        fixed_full,
        supplement;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :ggt_nearest,
        residual_keep_policy = :near_null_only,
    )
    legacy_ops = ordinary_cartesian_qiu_white_operators(
        fixed_legacy,
        supplement;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :ggt_nearest,
        residual_keep_policy = :near_null_only,
    )

    full_rep = basis_representation(full_ops)
    legacy_rep = basis_representation(legacy_ops)
    source_one_body = Matrix{Float64}(full_ops.one_body_hamiltonian)
    # Final working bases follow the repo orthonormal-basis contract, so the
    # occupied benchmark columns come from the ordinary Hermitian one-body
    # problem, not a generalized eigenproblem against `S = I + eps`.
    source_orbitals = eigen(Symmetric(source_one_body))
    occupied_order = sortperm(source_orbitals.values)
    source_occ = Matrix{Float64}(source_orbitals.vectors[:, occupied_order[1:1]])

    paths = _case_paths()
    mkpath(paths.dir)
    write_cartesian_basis_bundle_jld2(paths.source, full_ops; include_ham = true)
    write_cartesian_basis_bundle_jld2(paths.target, legacy_ops; include_ham = true)

    return (
        atom = "He",
        supplement_basis = "cc-pVTZ",
        supplement_lmax = 1,
        nside = 5,
        parent_count = 13,
        source_working_box = fixed_full.shell.working_box,
        target_working_box = fixed_legacy.shell.working_box,
        source_rep = full_rep,
        target_rep = legacy_rep,
        source_ops = full_ops,
        target_ops = legacy_ops,
        source_occ = source_occ,
        paths = paths,
    )
end

function _dense_overlap_stage(fixture)
    start_ns = time_ns()
    cross_overlap_matrix = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        fixture.target_rep,
        fixture.source_rep,
    )
    overlap_build_seconds = (time_ns() - start_ns) / 1.0e9
    return (
        transfer_matrix_size = size(cross_overlap_matrix),
        overlap_matrix_materialized = true,
        overlap_build_seconds = overlap_build_seconds,
        stage_timings = nothing,
    )
end

function _factorized_overlap_stage(fixture)
    start_ns = time_ns()
    projector_result = GaussletBases._cartesian_basis_projector_with_stage_timings(
        fixture.source_rep,
        fixture.target_rep,
    )
    overlap_build_seconds = (time_ns() - start_ns) / 1.0e9
    return (
        transfer_matrix_size = size(projector_result.projector.matrix),
        overlap_matrix_materialized = true,
        overlap_build_seconds = overlap_build_seconds,
        stage_timings = projector_result.stage_timings,
    )
end

function _dense_transfer_stage(fixture)
    # Final orthonormal-basis transfer uses only the exact cross overlap `S_BA`;
    # self-overlaps are diagnostic/debug data and are intentionally out of the
    # benchmark timing path.
    overlap_start_ns = time_ns()
    target_source_overlap = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        fixture.target_rep,
        fixture.source_rep,
    )
    overlap_build_seconds = (time_ns() - overlap_start_ns) / 1.0e9
    apply_start_ns = time_ns()
    coefficients = target_source_overlap * fixture.source_occ
    apply_seconds = (time_ns() - apply_start_ns) / 1.0e9
    return (
        coefficient_size = size(coefficients),
        overlap_matrix_materialized = true,
        reused_overlap_matrix = true,
        overlap_build_seconds = overlap_build_seconds,
        apply_seconds = apply_seconds,
        transfer_matrix_size = size(target_source_overlap),
        stage_timings = nothing,
    )
end

function _factorized_transfer_stage(fixture)
    projector_stage = _timed_stage(() ->
        GaussletBases._cartesian_basis_projector_with_stage_timings(
            fixture.source_rep,
            fixture.target_rep,
        )
    )
    projector_result = projector_stage.result
    materialized_stage =
        _timed_stage(() -> transfer_orbitals(fixture.source_occ, projector_result.projector))
    direct_stage = _timed_stage(() ->
        transfer_orbitals(
            fixture.source_occ,
            fixture.source_rep,
            fixture.target_rep;
            materialize_projector = false,
        )
    )
    materialized_result = materialized_stage.result
    direct_result = direct_stage.result
    return (
        coefficient_size = size(materialized_result.coefficients),
        overlap_matrix_materialized = true,
        reused_overlap_matrix = true,
        overlap_build_seconds = projector_stage.seconds,
        apply_seconds = materialized_stage.seconds,
        materialized_total_seconds = projector_stage.seconds + materialized_stage.seconds,
        no_projector_transfer_seconds = direct_stage.seconds,
        transfer_path = materialized_result.diagnostics.transfer_path,
        no_projector_transfer_path = direct_result.diagnostics.transfer_path,
        projector_materialized = materialized_result.projector !== nothing,
        no_projector_materialized = direct_result.projector !== nothing,
        coefficient_max_abs_residual = _coefficient_max_abs_residual(
            materialized_result.coefficients,
            direct_result.coefficients,
        ),
        transfer_matrix_size = size(projector_result.projector.matrix),
        projector_rss_before_bytes = projector_stage.rss_before_bytes,
        projector_rss_after_bytes = projector_stage.rss_after_bytes,
        materialized_apply_rss_before_bytes = materialized_stage.rss_before_bytes,
        materialized_apply_rss_after_bytes = materialized_stage.rss_after_bytes,
        no_projector_rss_before_bytes = direct_stage.rss_before_bytes,
        no_projector_rss_after_bytes = direct_stage.rss_after_bytes,
        stage_timings = projector_result.stage_timings,
    )
end

function _disk_factorized_transfer_stage(fixture)
    read_stage = _timed_stage(() -> (
        source_bundle = read_cartesian_basis_bundle(fixture.paths.source),
        target_bundle = read_cartesian_basis_bundle(fixture.paths.target),
    ))
    source_bundle = read_stage.result.source_bundle
    target_bundle = read_stage.result.target_bundle
    projector_stage = _timed_stage(() ->
        GaussletBases._cartesian_basis_projector_with_stage_timings(
            source_bundle.basis,
            target_bundle.basis,
        )
    )
    projector_result = projector_stage.result
    materialized_stage =
        _timed_stage(() -> transfer_orbitals(fixture.source_occ, projector_result.projector))
    bundle_direct_stage = _timed_stage(() ->
        transfer_orbitals(
            fixture.source_occ,
            source_bundle,
            target_bundle;
            materialize_projector = false,
        )
    )
    path_direct_stage = _timed_stage(() ->
        transfer_orbitals(
            fixture.source_occ,
            fixture.paths.source,
            fixture.paths.target;
            materialize_projector = false,
        )
    )
    materialized_result = materialized_stage.result
    bundle_direct_result = bundle_direct_stage.result
    path_direct_result = path_direct_stage.result
    return (
        coefficient_size = size(materialized_result.coefficients),
        overlap_matrix_materialized = true,
        reused_overlap_matrix = true,
        bundle_read_seconds = read_stage.seconds,
        overlap_build_seconds = projector_stage.seconds,
        apply_seconds = materialized_stage.seconds,
        materialized_total_seconds =
            read_stage.seconds + projector_stage.seconds + materialized_stage.seconds,
        no_projector_bundle_transfer_seconds = bundle_direct_stage.seconds,
        no_projector_path_transfer_seconds = path_direct_stage.seconds,
        transfer_path = materialized_result.diagnostics.transfer_path,
        no_projector_bundle_transfer_path = bundle_direct_result.diagnostics.transfer_path,
        no_projector_path_transfer_path = path_direct_result.diagnostics.transfer_path,
        projector_materialized = materialized_result.projector !== nothing,
        no_projector_bundle_materialized = bundle_direct_result.projector !== nothing,
        no_projector_path_materialized = path_direct_result.projector !== nothing,
        coefficient_max_abs_residual_bundle = _coefficient_max_abs_residual(
            materialized_result.coefficients,
            bundle_direct_result.coefficients,
        ),
        coefficient_max_abs_residual_path = _coefficient_max_abs_residual(
            materialized_result.coefficients,
            path_direct_result.coefficients,
        ),
        transfer_matrix_size = size(projector_result.projector.matrix),
        bundle_read_rss_before_bytes = read_stage.rss_before_bytes,
        bundle_read_rss_after_bytes = read_stage.rss_after_bytes,
        projector_rss_before_bytes = projector_stage.rss_before_bytes,
        projector_rss_after_bytes = projector_stage.rss_after_bytes,
        materialized_apply_rss_before_bytes = materialized_stage.rss_before_bytes,
        materialized_apply_rss_after_bytes = materialized_stage.rss_after_bytes,
        no_projector_bundle_rss_before_bytes = bundle_direct_stage.rss_before_bytes,
        no_projector_bundle_rss_after_bytes = bundle_direct_stage.rss_after_bytes,
        no_projector_path_rss_before_bytes = path_direct_stage.rss_before_bytes,
        no_projector_path_rss_after_bytes = path_direct_stage.rss_after_bytes,
        stage_timings = projector_result.stage_timings,
    )
end

function _mode_runner(mode::AbstractString)
    if mode == "dense_overlap"
        return _dense_overlap_stage
    elseif mode == "factorized_overlap"
        return _factorized_overlap_stage
    elseif mode == "dense_transfer"
        return _dense_transfer_stage
    elseif mode == "factorized_transfer"
        return _factorized_transfer_stage
    elseif mode == "disk_factorized_transfer"
        return _disk_factorized_transfer_stage
    end
    throw(
        ArgumentError(
            "unknown mode $(repr(mode)); expected dense_overlap, factorized_overlap, dense_transfer, factorized_transfer, or disk_factorized_transfer",
        ),
    )
end

function _measure_mode(fixture, mode::AbstractString; reps::Int)
    runner = _mode_runner(mode)

    runner(fixture)
    GC.gc()

    times = Float64[]
    summary = nothing
    for _ in 1:reps
        GC.gc()
        start_ns = time_ns()
        summary = runner(fixture)
        elapsed_seconds = (time_ns() - start_ns) / 1.0e9
        push!(times, elapsed_seconds)
    end

    return (
        case = "atomic_hybrid_he_ns5",
        mode = mode,
        reps = reps,
        wall_seconds = only(times),
        wall_seconds_samples = times,
        peak_rss_bytes = _rss_bytes(),
        atom = fixture.atom,
        supplement_basis = fixture.supplement_basis,
        supplement_lmax = fixture.supplement_lmax,
        nside = fixture.nside,
        parent_count = fixture.parent_count,
        source_working_box = fixture.source_working_box,
        target_working_box = fixture.target_working_box,
        source_final_dimension = fixture.source_rep.metadata.final_dimension,
        target_final_dimension = fixture.target_rep.metadata.final_dimension,
        source_raw_dimension = fixture.source_rep.metadata.parent_dimension,
        target_raw_dimension = fixture.target_rep.metadata.parent_dimension,
        source_gausslet_count = fixture.source_rep.metadata.route_metadata.gausslet_count,
        target_gausslet_count = fixture.target_rep.metadata.route_metadata.gausslet_count,
        source_residual_count = fixture.source_rep.metadata.route_metadata.residual_count,
        target_residual_count = fixture.target_rep.metadata.route_metadata.residual_count,
        source_occ_columns = size(fixture.source_occ, 2),
        summary = summary,
    )
end

function main(args)
    values = _parse_args(args)
    mode = get(values, "mode", "factorized_transfer")
    reps = parse(Int, get(values, "reps", "1"))
    reps > 0 || throw(ArgumentError("reps must be positive"))

    fixture = _he_ns5_fixture()
    println(_measure_mode(fixture, mode; reps = reps))
end

main(ARGS)
