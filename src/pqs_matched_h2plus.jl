"""
    PQSH2PlusRow

One result row from [`pqs_h2plus_comparison`](@ref). Public fields identify the
basis, dimension, captured parent-state norm, electronic and nuclear-repulsion
energies, total energy, and the contraction and independent-reference errors,
all in hartree where applicable.
"""
struct PQSH2PlusRow
    basis::Symbol
    dimension::Int
    capture::Float64
    electronic_energy_hartree::Float64
    nuclear_repulsion_energy_hartree::Float64
    total_energy_hartree::Float64
    contraction_error_hartree::Float64
    total_error_hartree::Float64
end

"""
    PQSH2PlusComparison

Matched full-parent, PQS, and White-Lindsey H2+ rows returned by
[`pqs_h2plus_comparison`](@ref), together with the supplied independent
reference energy and common parent/core/shell/slab dimension accounting.
"""
struct PQSH2PlusComparison
    rows::NTuple{3,PQSH2PlusRow}
    independent_reference_total_energy_hartree::Float64
    parent_axis_counts::NTuple{3,Int}
    parent_dimension::Int
    direct_core_columns::Int
    complete_shell_count::Int
    complete_shell_columns::Int
    slab_count::Int
    slab_columns::Int
end

_pqs_h2plus_axis_data(axis) = (
    axis.basis.reference_center_data,
    axis.basis.center_data,
    axis.basis.integral_weight_data,
    axis.basis.coefficient_matrix,
    axis.centers,
    axis.weights,
    axis.overlap,
    axis.kinetic,
    axis.position,
    axis.x2,
)

function _pqs_h2plus_parent_fingerprint(pgdg)
    bytes = vcat((reinterpret(UInt8, vec(array)) for axis in pgdg
        for array in _pqs_h2plus_axis_data(axis))...)
    return bytes2hex(sha256(bytes))
end

function _pqs_h2plus_mode_product!(out, storage, matrix, tensor, axis)
    order = axis == 1 ? (1, 2, 3) : axis == 2 ? (2, 1, 3) : (3, 1, 2)
    moved = PermutedDimsArray(tensor, order)
    product = reshape(storage, size(matrix, 1), :)
    mul!(product, matrix, reshape(moved, size(tensor, axis), :))
    shaped = reshape(product, size(matrix, 1), size(tensor, order[2]),
        size(tensor, order[3]))
    return permutedims!(out, shaped, invperm(order))
end

function _pqs_h2plus_product_apply!(out, matrices, vector, dimensions, storage)
    shape = reverse(dimensions)
    _pqs_h2plus_mode_product!(reshape(storage[1], shape), storage[3], matrices[3],
        reshape(vector, shape), 1)
    _pqs_h2plus_mode_product!(reshape(storage[2], shape), storage[3], matrices[2],
        reshape(storage[1], shape), 2)
    _pqs_h2plus_mode_product!(reshape(out, shape), storage[3], matrices[1],
        reshape(storage[2], shape), 3)
    return out
end

function _pqs_h2plus_inverse_sqrt(overlap)
    decomposition = eigen(Symmetric(overlap))
    minimum(decomposition.values) > 0 ||
        throw(ArgumentError("H2+ parent overlap is not positive definite"))
    root = decomposition.vectors * Diagonal(inv.(sqrt.(decomposition.values))) *
        transpose(decomposition.vectors)
    error = norm(root * overlap * root - I, Inf)
    error <= 1.0e-10 ||
        throw(ArgumentError("H2+ parent overlap orthogonalization failed"))
    return root
end

_pqs_h2plus_factor_term(factors, term) =
    factors isa AbstractArray{<:Real,3} ? @view(factors[term, :, :]) : factors[term]

function _pqs_h2plus_parent_solution(pgdg, expansion, nuclei)
    dimensions = ntuple(axis -> length(pgdg[axis].centers), 3)
    prod(dimensions) == 12789 ||
        throw(DimensionMismatch("matched H2+ parent dimension must be 12789"))
    overlaps = ntuple(axis -> pgdg[axis].overlap, 3)
    kinetic = ntuple(axis -> pgdg[axis].kinetic, 3)
    roots = ntuple(axis -> _pqs_h2plus_inverse_sqrt(overlaps[axis]), 3)
    factors = ntuple(center -> ntuple(axis ->
        _pqs_source_box_route_driver_centered_factor_terms(
            pgdg[axis], expansion, nuclei[center][axis]), 3), 2)
    for matrices in (overlaps, kinetic, roots), matrix in matrices
        all(isfinite, matrix) && norm(matrix - transpose(matrix), Inf) <= 1.0e-10 ||
            throw(ArgumentError("matched H2+ parent axis operator is invalid"))
    end
    for center in factors, axis in center, term in eachindex(expansion.coefficients)
        matrix = _pqs_h2plus_factor_term(axis, term)
        all(isfinite, matrix) && norm(matrix - transpose(matrix), Inf) <= 1.0e-10 ||
            throw(ArgumentError("matched H2+ parent nuclear factor is invalid"))
    end
    n = prod(dimensions)
    storage = ntuple(_ -> zeros(Float64, n), 3)
    orthogonal, raw, nuclear, product = ntuple(_ -> zeros(Float64, n), 4)
    product_apply!(out, matrices, vector) =
        _pqs_h2plus_product_apply!(out, matrices, vector, dimensions, storage)
    function raw_kinetic!(out, vector)
        product_apply!(out, (kinetic[1], overlaps[2], overlaps[3]), vector)
        for matrices in ((overlaps[1], kinetic[2], overlaps[3]),
            (overlaps[1], overlaps[2], kinetic[3]))
            product_apply!(product, matrices, vector)
            out .+= product
        end
        return out
    end
    function raw_nuclear!(out, vector, center)
        fill!(out, 0.0)
        for term in eachindex(expansion.coefficients)
            product_apply!(product, ntuple(axis ->
                _pqs_h2plus_factor_term(factors[center][axis], term), 3), vector)
            out .-= expansion.coefficients[term] .* product
        end
        return out
    end
    function apply_h!(out, vector)
        Base.mightalias(out, vector) &&
            throw(ArgumentError("matched H2+ parent action input and output alias"))
        product_apply!(orthogonal, roots, vector)
        raw_kinetic!(raw, orthogonal)
        for center in 1:2
            raw_nuclear!(nuclear, orthogonal, center)
            raw .+= nuclear
        end
        return product_apply!(out, roots, raw)
    end
    solved = _lanczos_ground_state_apply(apply_h!, prod(dimensions); tol = 1.0e-10)
    orbital = solved.vector
    reapplied = similar(orbital)
    apply_h!(reapplied, orbital)
    residual = norm(reapplied - solved.value .* orbital)
    residual <= 1.0e-9 ||
        throw(ArgumentError("matched H2+ parent eigen-residual exceeds 1.0e-9"))
    probe_a = normalize!(sin.(0.017 .* (1:prod(dimensions))) .+
        0.31 .* cos.(0.029 .* (1:prod(dimensions))))
    probe_b = normalize!(cos.(0.023 .* (1:prod(dimensions))) .-
        0.27 .* sin.(0.037 .* (1:prod(dimensions))))
    h_a, h_b = similar(probe_a), similar(probe_b)
    apply_h!(h_a, probe_a)
    apply_h!(h_b, probe_b)
    abs(dot(probe_a, h_b) - dot(h_a, probe_b)) <= 1.0e-10 ||
        throw(ArgumentError("matched H2+ parent H1 symmetry check failed"))
    return solved.value, orbital, overlaps, roots
end

function _pqs_h2plus_terminal_solution(base)
    products = cartesian_base_products(base)
    unit_nuclear = cartesian_base_unit_nuclear(base)
    H1 = copy(products.kinetic)
    for matrix in unit_nuclear
        H1 .+= matrix
    end
    all(isfinite, H1) ||
        throw(ArgumentError("matched H2+ terminal H1 contains nonfinite values"))
    norm(H1 - transpose(H1), Inf) <= 1.0e-10 ||
        throw(ArgumentError("matched H2+ terminal H1 symmetry check failed"))
    solved = eigen(Symmetric(H1))
    orbital = solved.vectors[:, 1]
    residual = norm(H1 * orbital - solved.values[1] .* orbital)
    residual <= 1.0e-9 ||
        throw(ArgumentError("matched H2+ terminal eigen-residual exceeds 1.0e-9"))
    return solved.values[1], orbital
end

function _pqs_h2plus_only(values, label)
    distinct = unique(values)
    length(distinct) == 1 ||
        throw(ArgumentError("matched H2+ region has inconsistent $(label)"))
    return only(distinct)
end

function _pqs_h2plus_region_signature(base, key)
    rows = base.terminal_due_diligence.terminal_rows
    selected = findall(row -> row.region_key === key, rows)
    isempty(selected) && throw(ArgumentError("matched H2+ region is empty"))
    blocks = base.terminal_basis.blocks[selected]
    index_bounds = ntuple(axis -> begin
        field = (:x, :y, :z)[axis]
        ranges = [getproperty(rows[i].index_ranges, field) for i in selected]
        minimum(first, ranges):maximum(last, ranges)
    end, 3)
    physical_bounds = ntuple(axis -> begin
        field = (:x, :y, :z)[axis]
        ranges = [getproperty(rows[i].physical_ranges, field) for i in selected]
        (minimum(first, ranges), maximum(last, ranges))
    end, 3)
    support = sort!(reduce(vcat, (block.support_indices for block in blocks)))
    return (
        _pqs_h2plus_only([rows[i].region_kind for i in selected], "kind"),
        _pqs_h2plus_only([rows[i].role for i in selected], "role"),
        _pqs_h2plus_only([rows[i].owner_contact_shared for i in selected], "owner"),
        _pqs_h2plus_only([rows[i].shell_index for i in selected], "shell index"),
        index_bounds, physical_bounds, support, sum(rows[i].retained_count for i in selected),
        _pqs_h2plus_only([rows[i].source_mode_shape for i in selected], "source shape"),
        _pqs_h2plus_only([(rows[i].slab_axis, rows[i].slab_side,
            rows[i].slab_thickness, rows[i].slab_stack_index,
            rows[i].slab_stack_count) for i in selected], "slab metadata"),
    )
end

function _pqs_h2plus_validate_route(base, nesting)
    input = base.input
    q = nesting === :pqs ? 5 : 3
    (input.kind, input.ns, input.q, input.core_spacing, input.s_factor,
        input.tail_spacing, input.xmax_parallel, input.xmax_transverse,
        input.source_span, input.coulomb_accuracy) ==
        (:h2, 5, q, 0.30, 1.0, 2.8, 11.0, 10.0, :ordinary, :high) ||
        throw(ArgumentError("matched H2+ resolved route inputs differ from the frozen contract"))
    length(base.coulomb_expansion.coefficients) == 135 ||
        throw(ArgumentError("matched H2+ requires the high135 Coulomb expansion"))
    basis = base.terminal_basis
    due = base.terminal_due_diligence
    basis.final_dimension == 1285 ||
        throw(DimensionMismatch("matched H2+ terminal dimension must be 1285"))
    rows = due.terminal_rows
    length(rows) == length(basis.blocks) ||
        throw(DimensionMismatch("matched H2+ terminal rows and blocks differ"))
    all(i -> rows[i].terminal_order == i &&
        rows[i].support_rows == length(basis.blocks[i].support_indices) &&
        rows[i].final_column_range == basis.blocks[i].column_range,
        eachindex(rows)) ||
        throw(ArgumentError("matched H2+ terminal row ordering differs from native blocks"))
    collect(Iterators.flatten(row.final_column_range for row in rows)) ==
        collect(1:1285) ||
        throw(ArgumentError("matched H2+ terminal columns do not partition the basis"))
    dimensions = due.dimensions
    (dimensions.direct_core_columns, dimensions.complete_shell_columns,
        dimensions.slab_columns) == (275, 960, 50) ||
        throw(DimensionMismatch("matched H2+ terminal column accounting differs"))
    length(unique(row.region_key for row in rows
        if row.region_kind === :complete_shell)) == 8 ||
        throw(ArgumentError("matched H2+ requires eight complete shells"))
    length(unique(row.region_key for row in rows
        if occursin("slab", String(row.region_kind)))) == 2 ||
        throw(ArgumentError("matched H2+ requires two slabs"))
    return nothing
end

function _pqs_h2plus_validate_common_parent(pqs, wl)
    p_pgdg, w_pgdg = _cartesian_base_pgdg(pqs), _cartesian_base_pgdg(wl)
    axis_counts = ntuple(axis -> length(p_pgdg[axis].centers), 3)
    axis_counts == (21, 21, 29) ||
        throw(DimensionMismatch("matched H2+ parent axes must be 21 x 21 x 29"))
    all(_pqs_h2plus_axis_data(left) == _pqs_h2plus_axis_data(right)
        for (left, right) in zip(p_pgdg, w_pgdg)) ||
        throw(ArgumentError("matched H2+ routes do not share one numerical parent"))
    _pqs_h2plus_parent_fingerprint(p_pgdg) ==
        _pqs_h2plus_parent_fingerprint(w_pgdg) ||
        throw(ArgumentError("matched H2+ parent fingerprints differ"))
    p_rows = pqs.terminal_due_diligence.terminal_rows
    w_rows = wl.terminal_due_diligence.terminal_rows
    p_keys = unique(row.region_key for row in p_rows)
    w_keys = unique(row.region_key for row in w_rows)
    p_keys == w_keys ||
        throw(ArgumentError("matched H2+ physical terminal region keys differ"))
    all(_pqs_h2plus_region_signature(pqs, key) ==
        _pqs_h2plus_region_signature(wl, key) for key in p_keys) ||
        throw(ArgumentError("matched H2+ physical terminal regions or supports differ"))
    for key in p_keys
        signature = _pqs_h2plus_region_signature(pqs, key)
        signature[1] === :complete_shell || continue
        shape = signature[9]
        shape[1] == shape[2] == 5 ||
            throw(ArgumentError("matched White-Lindsey shell outer shape differs from PQS"))
        inner = ntuple(axis -> shape[axis] - 2, 3)
        _pqs_h2plus_region_signature(wl, key)[8] == prod(shape) - prod(inner) ||
            throw(ArgumentError("matched White-Lindsey shell did not use axis-specific inner counts"))
    end
    return p_pgdg, axis_counts
end

function _pqs_h2plus_capture(base, parent_orbital, overlaps, roots, dimensions)
    blocks = base.terminal_basis.blocks
    supports = reduce(vcat, (block.support_indices for block in blocks))
    sort!(supports) == collect(1:prod(dimensions)) ||
        throw(ArgumentError("matched H2+ terminal supports do not partition the parent"))
    sqrt_factors = ntuple(axis -> overlaps[axis] * roots[axis], 3)
    n = prod(dimensions)
    storage = ntuple(_ -> zeros(Float64, n), 3)
    metric_orbital = zeros(Float64, n)
    _pqs_h2plus_product_apply!(metric_orbital, sqrt_factors, parent_orbital,
        dimensions, storage)
    parent_coefficients = zeros(Float64, n)
    capture = 0.0
    for block in blocks
        indices = block.support_indices
        source = @view metric_orbital[indices]
        amplitudes = isnothing(block.coefficients) ? source :
            transpose(block.coefficients) * source
        capture += sum(abs2, amplitudes)
        parent_coefficients[indices] .= isnothing(block.coefficients) ? amplitudes :
            block.coefficients * amplitudes
    end
    _pqs_h2plus_product_apply!(metric_orbital, sqrt_factors, parent_coefficients,
        dimensions, storage)
    residual = parent_orbital - metric_orbital
    _pqs_h2plus_product_apply!(metric_orbital, sqrt_factors, residual,
        dimensions, storage)
    orthogonality = 0.0
    for block in blocks
        indices = block.support_indices
        source = @view metric_orbital[indices]
        amplitudes = isnothing(block.coefficients) ? source :
            transpose(block.coefficients) * source
        orthogonality += sum(abs2, amplitudes)
    end
    parent_norm = dot(parent_orbital, parent_orbital)
    lost = sum(abs2, residual)
    abs(capture + lost - parent_norm) <= 1.0e-9 ||
        throw(ArgumentError("matched H2+ capture norm does not close"))
    sqrt(orthogonality) <= 1.0e-9 ||
        throw(ArgumentError("matched H2+ capture residual is not terminal-orthogonal"))
    return capture / parent_norm
end

"""
    pqs_h2plus_comparison(; independent_reference_total_energy_hartree)

Reconstruct the frozen R=2 H2+ full-parent, PQS, and matched White-Lindsey
one-body comparison. The external reference is used only for reported errors.
"""
function pqs_h2plus_comparison(;
    independent_reference_total_energy_hartree::Real,
)
    reference = Float64(independent_reference_total_energy_hartree)
    isfinite(reference) ||
        throw(ArgumentError("independent H2+ reference energy must be finite"))
    nuclei = [(0.0, 0.0, -1.0), (0.0, 0.0, 1.0)]
    system = (; atom_symbols = ["H", "H"], nuclear_charges = [1.0, 1.0],
        atom_locations = nuclei, nup = 1, ndn = 0)
    basis(nesting) = (; ns = 5, nesting, core_spacing = 0.30,
        xmax_parallel = 11.0, xmax_transverse = 10.0, s_factor = 1.0,
        tail_spacing = 2.8, source_span = :ordinary, coulomb_accuracy = :high)
    pqs = cartesian_base_working_basis(system; basis = basis(:pqs))
    wl = cartesian_base_working_basis(system; basis = basis(:wl))
    _pqs_h2plus_validate_route(pqs, :pqs)
    _pqs_h2plus_validate_route(wl, :wl)
    pgdg, axis_counts = _pqs_h2plus_validate_common_parent(pqs, wl)
    parent_energy, parent_orbital, overlaps, roots =
        _pqs_h2plus_parent_solution(pgdg, pqs.coulomb_expansion, nuclei)
    pqs_energy, _ = _pqs_h2plus_terminal_solution(pqs)
    wl_energy, _ = _pqs_h2plus_terminal_solution(wl)
    parent_total = parent_energy + 0.5
    pqs_total = pqs_energy + 0.5
    wl_total = wl_energy + 0.5
    min(pqs_total, wl_total) + 1.0e-10 >= parent_total ||
        throw(ArgumentError("matched H2+ terminal energy is below the full parent"))
    dimensions = axis_counts
    pqs_capture = _pqs_h2plus_capture(pqs, parent_orbital, overlaps, roots, dimensions)
    wl_capture = _pqs_h2plus_capture(wl, parent_orbital, overlaps, roots, dimensions)
    rows = (
        PQSH2PlusRow(:parent, 12789, 1.0, parent_energy, 0.5,
            parent_total, 0.0, parent_total - reference),
        PQSH2PlusRow(:pqs, 1285, pqs_capture, pqs_energy, 0.5,
            pqs_total, pqs_total - parent_total, pqs_total - reference),
        PQSH2PlusRow(:white_lindsey, 1285, wl_capture, wl_energy, 0.5,
            wl_total, wl_total - parent_total, wl_total - reference),
    )
    return PQSH2PlusComparison(rows, reference, axis_counts, 12789,
        275, 8, 960, 2, 50)
end
