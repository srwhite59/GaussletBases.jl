using GaussletBases
using JLD2
using LinearAlgebra
using Test

const H_LOWEST = -0.49877574806444014
const H2_LOWEST = -0.79460371733658908
const H2_SELF_COULOMB = 0.4569012290840094
const ATOL = 1.0e-10
const PROVENANCE_KEYS = (
    :provenance_version, :producer, :route, :q, :core_spacing,
    :reference_spacing, :tail_spacing, :parent_axis_family, :parent_axis_counts,
    :s_factor, :mapping_kind, :mapping_d, :mapping_s_factor,
    :mapping_s_standard, :mapping_s_effective,
    :radius, :xmax_parallel, :xmax_transverse,
    :atom_symbols, :nuclear_charges, :atom_locations, :nup, :ndn,
    :final_dimension,
)

const H_BASIS = (;
    q = 5,
    core_spacing = 0.5,
    radius = 4.0,
    reference_spacing = 1.0,
    d = 0.5,
)

const H2_BASIS = (;
    ns = 5,
    core_spacing = 0.5,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
)

const H_ACCURACY_BASIS = (;
    ns = 3,
    core_spacing = 0.5,
    radius = 2.0,
    reference_spacing = 1.0,
)

h_system() = (;
    atom_symbols = ["H"],
    nuclear_charges = [1.0],
    atom_locations = [(0.0, 0.0, 0.0)],
    nup = 1,
    ndn = 0,
)

h2_system() = (;
    atom_symbols = ["H", "H"],
    nuclear_charges = [1.0, 1.0],
    atom_locations = [(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)],
    nup = 1,
    ndn = 1,
)

function lowest_one_body(ham)
    H = one_body_hamiltonian(ham)
    eig = eigen(Symmetric((H + transpose(H)) ./ 2))
    return minimum(eig.values), eig.vectors[:, argmin(eig.values)], H
end

function self_coulomb(V, orbital)
    density = orbital * transpose(orbital)
    rho = 0.5 .* (density .+ transpose(density))
    occupations = vec(diag(rho))
    sym = 0.5 .* (V .+ transpose(V))
    return 2.0 * dot(occupations, sym * occupations) -
        dot(vec(rho), vec(sym .* rho))
end

function check_finite_symmetric(matrix)
    @test all(isfinite, matrix)
    @test norm(matrix - transpose(matrix), Inf) <= ATOL
end

function check_sector_parity(left_base, left, right_base, right)
    for axis in (:x, :y, :z)
        lhs = GaussletBases._nested_axis_pgdg(left_base.parent.parent_axis_bundle_object, axis)
        rhs = GaussletBases._nested_axis_pgdg(right_base.parent.parent_axis_bundle_object, axis)
        @test lhs.centers == rhs.centers
        @test lhs.weights == rhs.weights
    end
    block_signature(base) = [(block.support_indices, block.support_states,
        block.column_range, block.coefficients) for block in base.terminal_basis.blocks]
    @test block_signature(left_base) == block_signature(right_base)
    @test left.kinetic == right.kinetic
    @test left.nuclear_attraction_unit_by_center == right.nuclear_attraction_unit_by_center
    @test one_body_hamiltonian(left) == one_body_hamiltonian(right)
    @test left.electron_electron_ida == right.electron_electron_ida
    @test left.nuclear_repulsion == right.nuclear_repulsion
end
function check_provenance_keys(file)
    for key in PROVENANCE_KEYS
        @test haskey(file, "producer_provenance/$(key)")
    end
end

function check_coulomb_summary(file, policy, terms)
    @test file["coulomb_expansion/policy"] === policy
    @test file["coulomb_expansion/doacc"] === (policy === :high)
    @test file["coulomb_expansion/term_count"] == terms
    expected = policy === :high ? (1.0, 0.16, 0.01, 135.0) :
        (0.6, 0.5, 0.03, 27.0)
    @test Tuple(file["coulomb_expansion/$(key)"] for key in (:del, :s, :c, :maxu)) ==
        expected
end

function accuracy_case(basis)
    base = GaussletBases.cartesian_base_working_basis(h_system(); basis)
    ham = GaussletBases.cartesian_base_hamiltonian_assembly(base)
    pgdg = Tuple(GaussletBases._nested_axis_pgdg(
        base.parent.parent_axis_bundle_object, axis) for axis in (:x, :y, :z))
    return (; base, ham, pgdg)
end

function due_summary(case)
    due = case.base.terminal_due_diligence
    rows = isempty(due.terminal_rows) ? case.base.terminal_inventory.rows :
        due.terminal_rows
    retained(row) = hasproperty(row, :retained_count) ? row.retained_count :
        row.final_cols
    return (; bounds = due.geometry.parent_physical_bounds,
        axes = due.geometry.parent_axis_counts,
        radius = due.geometry.radius,
        padding = due.geometry.xmax_transverse,
        final_dimension = due.dimensions.base_final_dimension,
        retained = sum((retained(row) for row in rows); init = 0),
        topology = unique([row.region_kind for row in rows]),
        warnings = due.warnings)
end

@testset "public Cartesian base Hamiltonian" begin
    omitted_timed = @timed accuracy_case(H_ACCURACY_BASIS)
    compact_timed = @timed accuracy_case(
        merge(H_ACCURACY_BASIS, (; coulomb_accuracy = :compact)))
    omitted, compact = omitted_timed.value, compact_timed.value
    @test omitted.base.input.coulomb_accuracy === :compact
    @test length(omitted.base.coulomb_expansion) == 45
    @test omitted.ham.kinetic == compact.ham.kinetic
    @test omitted.ham.nuclear_attraction_unit_by_center ==
        compact.ham.nuclear_attraction_unit_by_center
    @test omitted.ham.electron_electron_ida == compact.ham.electron_electron_ida

    compact_factors = Tuple(axis.gaussian_factor_terms for axis in compact.pgdg)
    for invalid in (NaN, Inf, -Inf)
        destination = fill(0.25, compact.base.terminal_basis.final_dimension,
            compact.base.terminal_basis.final_dimension)
        unchanged = copy(destination)
        coefficients = copy(compact.base.coulomb_expansion.coefficients)
        coefficients[1] = invalid
        @test_throws ArgumentError begin
            GaussletBases.CartesianFinalBasisRealization._accumulate_terminal_gaussian_sum!(
                destination, compact.base.terminal_basis, coefficients,
                compact_factors[1], compact_factors[2], compact_factors[3])
        end
        @test destination == unchanged
    end

    high_timed = @timed accuracy_case(
        merge(H_ACCURACY_BASIS, (; coulomb_accuracy = :high)))
    high = high_timed.value
    @test high.base.input.coulomb_accuracy === :high
    @test length(high.base.coulomb_expansion) == 135
    @test all(axis -> axis.exponents == high.base.coulomb_expansion.exponents,
        high.pgdg)
    check_finite_symmetric(high.ham.kinetic)
    check_finite_symmetric(high.ham.electron_electron_ida)
    foreach(check_finite_symmetric, high.ham.nuclear_attraction_unit_by_center)

    wl_basis = merge(H_ACCURACY_BASIS,
        (; ns = 5, nesting = :wl, coulomb_accuracy = :high))
    wl = accuracy_case(wl_basis)
    @test wl.base.input.coulomb_accuracy === :high
    @test all(axis -> axis.exponents == wl.base.coulomb_expansion.exponents,
        wl.pgdg)
    check_finite_symmetric(wl.ham.electron_electron_ida)

    mktempdir() do dir
        high_path = joinpath(dir, "h_high_cartesian_ida.jld2")
        GaussletBases._cartesian_base_write_hamiltonian(
            high_path, high.ham, high.base)
        jldopen(high_path, "r") do file
            check_coulomb_summary(file, :high, 135)
            @test !haskey(file, "recipe_provenance/coulomb_accuracy")
        end
    end
    println("coulomb_accuracy_timing compact_s=", compact_timed.time,
        " compact_bytes=", compact_timed.bytes,
        " high_s=", high_timed.time, " high_bytes=", high_timed.bytes)
    println("coulomb_accuracy_due_compact=", due_summary(compact))
    println("coulomb_accuracy_due_high=", due_summary(high))
    println("coulomb_accuracy_due_wl=", due_summary(wl))

    h = cartesian_base_hamiltonian(h_system(); basis = H_BASIS)
    h_lowest, _, _ = lowest_one_body(h)
    @test h isa CartesianIDAHamiltonian{Float64}
    @test h_lowest ≈ H_LOWEST atol = ATOL
    check_finite_symmetric(h.kinetic)
    check_finite_symmetric(h.electron_electron_ida)
    foreach(check_finite_symmetric, h.nuclear_attraction_unit_by_center)

    he = merge(h_system(), (;
        atom_symbols = ["He"], nuclear_charges = [2.0], nup = 1, ndn = 1))
    he_plus = merge(he, (; nup = 1, ndn = 0))
    he_base = GaussletBases.cartesian_base_working_basis(he; basis = H_ACCURACY_BASIS)
    he_plus_base = GaussletBases.cartesian_base_working_basis(he_plus; basis = H_ACCURACY_BASIS)
    he_ham = GaussletBases.cartesian_base_hamiltonian_assembly(he_base)
    he_plus_ham = GaussletBases.cartesian_base_hamiltonian_assembly(he_plus_base)
    check_sector_parity(he_base, he_ham, he_plus_base, he_plus_ham)
    h2_base = GaussletBases.cartesian_base_working_basis(h2_system(); basis = H2_BASIS)
    h2 = GaussletBases.cartesian_base_hamiltonian_assembly(h2_base)
    h2_lowest, orbital, h2_one_body = lowest_one_body(h2)
    @test h2 isa CartesianIDAHamiltonian{Float64}
    @test size(h2.kinetic) == (487, 487)
    @test h2_lowest ≈ H2_LOWEST atol = ATOL
    @test self_coulomb(h2.electron_electron_ida, orbital) ≈ H2_SELF_COULOMB atol = ATOL
    check_finite_symmetric(h2.kinetic)
    check_finite_symmetric(h2.electron_electron_ida)
    foreach(check_finite_symmetric, h2.nuclear_attraction_unit_by_center)

    he2 = merge(h2_system(), (;
        atom_symbols = ["He", "He"], nuclear_charges = [2.0, 2.0], nup = 2, ndn = 2))
    he2_2plus = merge(he2, (; nup = 1, ndn = 1))
    sector_basis = merge(H2_BASIS, (; ns = 3))
    he2_base = GaussletBases.cartesian_base_working_basis(he2; basis = sector_basis)
    he2_2plus_base = GaussletBases.cartesian_base_working_basis(he2_2plus; basis = sector_basis)
    he2_ham = GaussletBases.cartesian_base_hamiltonian_assembly(he2_base)
    he2_2plus_ham = GaussletBases.cartesian_base_hamiltonian_assembly(he2_2plus_base)
    check_sector_parity(he2_base, he2_ham, he2_2plus_base, he2_2plus_ham)
    wl_h2_base = GaussletBases.cartesian_base_working_basis(h2_system(); basis = merge(H2_BASIS, (; nesting = :wl)))
    wl_h2 = GaussletBases.cartesian_base_hamiltonian_assembly(wl_h2_base)
    @test all(axis -> let p = GaussletBases._nested_axis_pgdg(h2_base.parent.parent_axis_bundle_object, axis),
            w = GaussletBases._nested_axis_pgdg(wl_h2_base.parent.parent_axis_bundle_object, axis)
        p.centers == w.centers && p.weights == w.weights end, (:x, :y, :z))
    pqs_shells = filter(row -> row.region_kind === :complete_shell, h2_base.terminal_due_diligence.terminal_rows)
    wl_shells = filter(row -> row.region_kind === :complete_shell, wl_h2_base.terminal_due_diligence.terminal_rows)
    for shell in pqs_shells
        children = filter(row -> row.region_key === shell.region_key, wl_shells)
        @test sum(row.retained_count for row in children) == shell.retained_count ==
            prod(shell.source_mode_shape) - prod(dim - 2 for dim in shell.source_mode_shape)
        @test sum(row.support_rows for row in children) == shell.support_rows
        for child in children
            free_axes = findall(>(1), child.outer_shape)
            @test child.retained_count == prod((shell.source_mode_shape[axis] - 2 for axis in free_axes); init = 1)
        end
    end
    @test wl_h2_base.terminal_basis.final_dimension == h2_base.terminal_basis.final_dimension
    @test collect(Iterators.flatten(block.column_range for block in wl_h2_base.terminal_basis.blocks)) == collect(1:wl_h2_base.terminal_basis.final_dimension)
    overlaps = Tuple(GaussletBases._nested_axis_pgdg(wl_h2_base.parent.parent_axis_bundle_object, axis).overlap for axis in (:x, :y, :z))
    wl_gram_error = maximum(block -> isnothing(block.coefficients) ? 0.0 :
        norm(transpose(block.coefficients) * GaussletBases.CartesianFinalBasisRealization._support_action(
            block.support_states, block.support_states, block.coefficients, overlaps) - I, Inf),
        wl_h2_base.terminal_basis.blocks)
    @test wl_gram_error <= ATOL
    foreach(check_finite_symmetric, (wl_h2.kinetic, wl_h2.electron_electron_ida, wl_h2.nuclear_attraction_unit_by_center...))

    mktempdir() do dir
        h_path = joinpath(dir, "h_cartesian_ida.jld2")
        written_h = cartesian_base_hamiltonian(h_system(); basis = H_BASIS, hamfile = h_path)
        jldopen(h_path, "r") do file
            check_provenance_keys(file)
            check_coulomb_summary(file, :compact, 45)
            @test !haskey(file, "recipe_provenance/coulomb_accuracy")
            @test file["producer_provenance/reference_spacing"] == 1.0
            @test file["producer_provenance/mapping_d"] == 0.5
            @test file["producer_provenance/s_factor"] == 1.0
            @test file["producer_provenance/mapping_s_factor"] == 1.0
            @test file["producer_provenance/mapping_s_standard"] == sqrt(0.5)
            @test file["producer_provenance/mapping_s_effective"] == sqrt(0.5)
            @test file["producer_provenance/mapping_kind"] === :white_lindsey_atomic_mapping
            @test file["producer_provenance/route"] === :one_center_pqs_base
            @test file["producer_provenance/final_dimension"] == size(written_h.kinetic, 1)
        end
        path = joinpath(dir, "h2_cartesian_ida.jld2")
        written = cartesian_base_hamiltonian(h2_system(); basis = H2_BASIS, hamfile = path)
        readback = read_cartesian_ida_hamiltonian(path)
        @test norm(one_body_hamiltonian(written) - one_body_hamiltonian(readback), Inf) == 0.0
        @test norm(h2_one_body - one_body_hamiltonian(readback), Inf) == 0.0
        jldopen(path, "r") do file
            check_provenance_keys(file)
            check_coulomb_summary(file, :compact, 45)
            @test file["producer_provenance/mapping_kind"] === :multicenter_pqs_mapping
            @test file["producer_provenance/mapping_d"] === nothing
            @test file["producer_provenance/route"] === :z_axis_diatomic_pqs_base
        end
        charged_path = joinpath(dir, "he2_2plus_cartesian_ida.jld2")
        GaussletBases._cartesian_base_write_hamiltonian(
            charged_path, he2_2plus_ham, he2_2plus_base)
        charged_readback = read_cartesian_ida_hamiltonian(charged_path)
        @test (charged_readback.nup, charged_readback.ndn) == (1, 1)
        @test charged_readback.kinetic == he2_2plus_ham.kinetic
        @test charged_readback.nuclear_attraction_unit_by_center ==
            he2_2plus_ham.nuclear_attraction_unit_by_center
        @test charged_readback.electron_electron_ida == he2_2plus_ham.electron_electron_ida
    end

    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(h_system(), (; extra = true)); basis = H_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h_system(); basis = merge(H_BASIS, (; xmax_parallel = 6.0)))
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h_system(); basis = merge(H_BASIS, (; s_factor = 0.0)))
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h_system(); basis = merge(H_BASIS, (; coulomb_accuracy = :medium)))
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h2_system(); basis = merge(H2_BASIS, (; radius = 4.0)))
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h2_system(); basis = merge(H2_BASIS, (; d = 0.3)))
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h_system(); basis = merge(H_BASIS, (; parent_axis_family = :G8)))
    @test_throws ArgumentError cartesian_base_hamiltonian(
        (; atom_symbols = ("H",), nuclear_charges = [1.0],
            atom_locations = [(0.0, 0.0, 0.0)], nup = 1, ndn = 0);
        basis = H_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(h_system(), (; atom_locations = [[0.0, 0.0, 0.0]])); basis = H_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(h2_system(), (; atom_locations = [(-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)]));
        basis = H2_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(h2_system(), (; atom_locations = [(1.0, 0.0, -2.0), (1.0, 0.0, 2.0)]));
        basis = H2_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(h2_system(), (; atom_locations = [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]));
        basis = H2_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(h2_system(), (; atom_locations = [(0.0, 0.0, -2.0), (1.0, 0.0, 2.0)]));
        basis = H2_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h2_system(); basis = H2_BASIS, hamfile = "")
    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(he, (; nup = 0, ndn = 0)); basis = H_ACCURACY_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(he, (; nup = -1, ndn = 1)); basis = H_ACCURACY_BASIS)
end

@testset "Scientific collinear q" begin
    G = GaussletBases; C = G.CartesianFinalBasisRealization
    z = [-1.8, 0., 1.8]; expansion = coulomb_gaussian_expansion(doacc=false)
    frozen_L = ([5,5,5,5,4,4,4], [9,7,6,5,5,5,5,5], [9,8,7,7,7,7,6,6])
    for q in 4:6
        w = cartesian_collinear_working_basis(z, ones(3); q, expansion)
        pg = ntuple(a -> G._nested_axis_pgdg(w.parent_axis_bundles, (:x,:y,:z)[a]), 3)
        dims = length.(getproperty.(pg, :weights)); S = getproperty.(pg, :overlap)
        @test dims == ((19,19,29), (21,21,35), (23,23,41))[q-3]
        @test w.terminal_basis.final_dimension == (797,1371,2569)[q-3]
        counts = zeros(Int, prod(dims)); shell_counts = Int[]
        for b in w.terminal_basis.blocks
            counts[b.support_indices] .+= 1
            Q = isnothing(b.coefficients) ? Matrix{Float64}(I, length(b.support_states), length(b.support_states)) : b.coefficients
            metric = C._prf_terminal_overlap(w.terminal_basis, b.support_states, Q, S)
            metric[b.column_range, :] -= I
            @test maximum(abs, metric) <= 1e-10
            weights = Q' * [prod(pg[a].weights[s[a]] for a in 1:3) for s in b.support_states]
            @test all(x -> isfinite(x) && x > 1e-14, weights)
            startswith(string(b.unit_key), "shell_") && push!(shell_counts, size(Q, 2))
        end
        @test all(==(1), counts)
        selected = Int[]; retention = G._nested_resolve_complete_shell_retention(q)
        function shell!(outer, inner)
            plan = G._nested_diatomic_source_box_dimension_plan((; nuclei=[(0.,0.,p) for p in z]),
                w.parent_axis_bundles, outer, inner, retention; bond_axis=:z, nside=q,
                selected_q=q, shared_shell_angular_resolution_scale=1.4)
            push!(selected, plan.source_mode_dims[3])
            @test all((q,q,last(selected))[a] <= length(outer[a]) for a in 1:3)
        end
        G.CartesianShellification._collinear_terminal_geometry((a,b)->nothing,
            shell!, p->nothing, getproperty.(pg, :centers), z, isodd(q) ? q : q+1, q)
        @test selected == frozen_L[q-3]
        @test shell_counts == [q*q*L - (q-2)^2*(L-2) for L in selected]
        if q == 5
            fixed = cartesian_collinear_working_basis(z, ones(3); q=4,
                core_spacing=.3, transverse_spacing=.3, expansion)
            @test all(G._nested_axis_pgdg(fixed.parent_axis_bundles, a).centers == pg[i].centers
                for (i,a) in enumerate((:x,:y,:z)))
        end
    end
    for q in (true, 2, 4.0, Inf, NaN, big(typemax(Int))+1)
        @test_throws ArgumentError cartesian_collinear_working_basis(z, ones(3); q, expansion)
    end
    for kw in ((; core_spacing=.3), (; transverse_spacing=.3), (; core_side=7),
        (; angular_reference_count=4), (; angular_resolution_scale=1.5),
        (; padding_parallel=0), (; tail_spacing=Inf), (; outer_face_count=0))
        @test_throws ArgumentError cartesian_collinear_working_basis(z, ones(3); q=5, expansion, kw...)
    end
    @test_throws ArgumentError cartesian_collinear_working_basis(z, [1.,2.,1.]; q=5, expansion)
end

@testset "Finite collinear PQS" begin
    GB = GaussletBases; C = GB.CartesianFinalBasisRealization
    expansion = coulomb_gaussian_expansion(doacc = false)
    controls = (; core_spacing = .6, transverse_spacing = .6, padding_parallel = 3.,
        padding_transverse = 3., core_side = 3, angular_reference_count = 5,
        outer_face_count = 3, tail_spacing = 2.8, angular_resolution_scale = 1.4, expansion)
    build(z, Z; kw...) = cartesian_collinear_working_basis(z, Z; merge(controls, (; kw...))...)
    function dense_map(basis, n)
        Q = zeros(n, basis.final_dimension)
        for b in basis.blocks
            Q[b.support_indices, b.column_range] = isnothing(b.coefficients) ?
                Matrix{Float64}(I, length(b.support_indices), length(b.column_range)) : b.coefficients
        end
        return Q
    end
    function action(v, A, B, D)
        nx, ny, nz = size(A, 1), size(B, 1), size(D, 1)
        a = reshape(D * reshape(v, nz, :), nz, ny, nx)
        b = reshape(B * reshape(permutedims(a, (2, 1, 3)), ny, :), ny, nz, nx)
        c = reshape(A * reshape(permutedims(b, (3, 1, 2)), nx, :), nx, ny, nz)
        return vec(permutedims(c, (3, 2, 1)))
    end
    function oracle(Q, X, pgdg, z, Z)
        A = Q * X; Y = zeros(size(A))
        S = ntuple(i -> pgdg[i].overlap, 3); T = ntuple(i -> pgdg[i].kinetic, 3)
        for j in axes(X, 2), a in 1:3
            Y[:, j] += action(A[:, j], ntuple(i -> i == a ? T[i] : S[i], 3)...)
        end
        for (p, charge) in zip(z, Z)
            f = ntuple(i -> C._terminal_factor_terms(
                GB._pqs_source_box_route_driver_centered_factor_terms(pgdg[i], expansion, i == 3 ? p : 0.)), 3)
            for j in axes(X, 2), k in eachindex(expansion.coefficients)
                Y[:, j] -= charge * expansion.coefficients[k] * action(A[:, j], ntuple(i -> view(f[i], k, :, :), 3)...)
            end
        end
        w = Q' * kron(pgdg[1].weights, pgdg[2].weights, pgdg[3].weights)
        @test all(w .> 1e-14)
        d = Q * (abs2.(X) ./ w); J = zeros(size(d))
        f = ntuple(i -> C._terminal_factor_terms(pgdg[i].pair_factor_terms_raw), 3)
        for j in axes(X, 2), k in eachindex(expansion.coefficients)
            J[:, j] += expansion.coefficients[k] * action(d[:, j], ntuple(i -> view(f[i], k, :, :), 3)...)
        end
        return X' * X, A' * Y, vec(sum(d .* J; dims = 1))
    end
    frozen = [
        ([.9999318254176067, .9997143222148448, .999713629868302],
         [-1.0801855049207716, -.6346238526272604, -.6346287143129328],
         [.35703227688764216, .2922260619143036, .29222484997191683]),
        ([.9997552845017766, .9982742935176115, .9982742028183909],
         [-1.2389311176399211, -.40792441475173197, -.407925743648283],
         [.7138968614711348, .5824489517158751, .5824489088945304]),
        ([.9988312030379334, .9911419746798421, .9911418012690119],
         [-1.0333535251547958, -.6019525037733703, -.6019518302045384],
         [.35666597897711483, .2888277307699443, .28882756100434104]),
        ([.9992927651276439, .9977624440489838, .9977623957328754],
         [-.9527395973938828, -.2826529501398075, -.2826553910092518],
         [.7123320035040893, .580175574336601, .5801755502698565])]
    limits = [(1.1e-6, 4.4e-6, 9.3e-6), (1.7e-7, 3.0e-6, 6.3e-8),
        (8.9e-6, 1.1e-5, 6.1e-5), (1.2e-7, 3.4e-6, 4.7e-8)]
    fixtures = [([-2.4, 0., 2.4], ones(3), .45, 6., 9),
        ([-2.4, -1.2, 2.4], ones(3), .6, 3., 7),
        ([-2.4, 0., 2.4], ones(3), .6, 3., 3),
        ([-3.6, -1.2, 1.2, 3.6], ones(4), .6, 3., 3),
        ([-1.2, 0., 1.2], ones(3), .6, 3., 3),
        ([-2.4, 0., 2.4], [1., 2., 1.], .6, 3., 3)]
    for (fixture, (z, Z, spacing, padding, outer)) in enumerate(fixtures)
        w = build(z, Z; transverse_spacing = spacing, padding_transverse = padding,
            outer_face_count = outer, core_side = UInt(3),
            angular_reference_count = Int32(fixture <= 2 ? 5 : 3),
            angular_resolution_scale = fixture == 3 ? Float32(1.4) : 1.4)
        basis, bundles = w.terminal_basis, w.parent_axis_bundles
        pgdg = ntuple(i -> GB._nested_axis_pgdg(bundles, (:x, :y, :z)[i]), 3)
        dims = ntuple(i -> length(pgdg[i].weights), 3); n = prod(dims)
        counts = zeros(Int, n)
        for b in basis.blocks
            counts[b.support_indices] .+= 1
            @test issorted(b.support_indices)
        end
        @test all(==(1), counts)
        Q = dense_map(basis, n)
        @test maximum(abs, Q' * kron(pgdg[1].overlap, pgdg[2].overlap, pgdg[3].overlap) * Q - I) <= 1e-10
        ops = cartesian_collinear_operators(w, z, Z; expansion)
        @test keys(ops) == (:one_body, :electron_electron_ida, :nuclear_repulsion)
        @test all(isfinite, ops.one_body) && all(isfinite, ops.electron_electron_ida)
        @test maximum(abs, ops.one_body - ops.one_body') <= 1e-10
        @test maximum(abs, ops.electron_electron_ida - ops.electron_electron_ida') <= 1e-10
        @test ops.nuclear_repulsion ≈ sum(Z[i]*Z[j]/(z[j]-z[i]) for i in eachindex(z) for j in i+1:length(z)) atol=1e-12 rtol=0
        states = [GB._cartesian_unflat_index(i, dims) for i in 1:n]
        parent = C.CartesianTerminalBasisRealization([C.CartesianTerminalBasisBlock(:parent,
            collect(1:n), states, nothing, 1:n)], n, 0.)
        if fixture <= 2
            @test (dims, basis.final_dimension) == (fixture == 1 ? ((13,13,17),1265) : ((9,9,15),619))
            @test Base.summarysize((w, ops)) < 512*1024^2
        end
        for (a, alpha) in enumerate((.1, .4))
            probes = CartesianGaussianShellSupplementRepresentation3D(:test,
                [CartesianGaussianShellOrbitalRepresentation3D(string(p), p, (0.,0.,0.), [alpha], [1.],
                 :axiswise_normalized_cartesian_gaussian) for p in ((0,0,0),(1,0,0),(0,1,0))], (;))
            X = gto_overlap_matrix(w, probes)
            Xp = C._terminal_residual_mixed_overlap(parent, bundles, probes)
            @test maximum(abs, X-Q'*Xp) <= 1e-12
            @test gto_overlap_matrix(w, probes; block_indices = [1,3]) == X[[1,3],:]
            packet = ExternalGTOOrbitalPacket(probes, Matrix{Float64}(I,3,3),
                ExternalGTOOrbitalSpinBlock(:restricted, Matrix{Float64}(I,3,3), ones(3)))
            @test import_external_gto_orbitals(w, packet).alpha.imported_coefficients == X
            gram, H, selfs = oracle(Q, X, pgdg, z, Z)
            @test maximum(abs, X'*ops.one_body*X-H) <= 1e-10
            @test maximum(abs, vec(sum(abs2.(X).*(ops.electron_electron_ida*abs2.(X));dims=1))-selfs) <= 1e-10
            if fixture <= 2
                index = 2*(fixture-1)+a; target = frozen[index]
                @test maximum(abs, diag(gram)-target[1]) <= 1e-8
                @test maximum(abs, H-Diagonal(target[2])) <= 1e-8
                @test maximum(abs, selfs-target[3]) <= 1e-8
                # Direct completion is an oracle only; unchanged inner columns are reused.
                inner = filter(b -> !startswith(string(b.unit_key), "outer_"), basis.blocks)
                outer_rows = sort(vcat([b.support_indices for b in basis.blocks if startswith(string(b.unit_key), "outer_")]...))
                inner_n = last(last(inner).column_range)
                @test inner_n == (fixture == 1 ? 617 : 423)
                Qd = [Q[:,1:inner_n] Matrix{Float64}(I,n,n)[:,outer_rows]]
                gd, hd, sd = oracle(Qd, Qd'*Xp, pgdg, z, Z)
                @test maximum(abs, diag(gd)-diag(gram)) <= limits[index][1]
                @test maximum(abs, hd-H) <= limits[index][2]
                @test maximum(abs, sd-selfs) <= limits[index][3]
            end
        end
        @test_throws ArgumentError cartesian_collinear_operators(w, z, Z; expansion = coulomb_gaussian_expansion())
    end
    # Three surviving groups: test exact direct gaps and realize the exterior slabs.
    grid_axes = (collect(-3.:3.), collect(-3.:3.), collect(-12.:12.))
    dims = (7,7,25); coverage = zeros(Int, dims); gap_count = Ref(0); exterior = Ref(0)
    bs = [build_basis(MappedUniformBasisSpec(:G10; count=n, mapping=IdentityMapping())) for n in (7,25)]
    bundle(b) = GB._mapped_ordinary_gausslet_1d_bundle(b; exponents=expansion.exponents,
        center=0., backend=:pgdg_localized_experimental)
    bx, bz = bundle(bs[1]), bundle(bs[2]); bundles = GB._CartesianNestedAxisBundles3D(bx,bx,bz)
    function mark(box, excluded)
        ss = [(x,y,z) for x in box[1] for y in box[2] for z in box[3]
            if !any(all((x,y,z)[a] in b[a] for a in 1:3) for b in excluded)]
        for s in ss; coverage[s...] += 1; end
        return length(ss)
    end
    direct(box, excluded) = (n=mark(box,excluded); isempty(excluded) || (gap_count[]+=n))
    shell(outer, inner) = mark(outer,[inner])
    function slab(piece)
        exterior[] += mark(piece.box, [])
        indices, states, coefficients = C._terminal_compact_thin_slab_block(
            [GB.CartesianCPB.cpb(piece.box)],piece.metadata,bundles)
        S = ntuple(a -> GB._nested_axis_pgdg(bundles,(:x,:y,:z)[a]).overlap,3)
        @test size(coefficients) == (49,9)
        @test maximum(abs, coefficients'*C._support_action(states,states,coefficients,S)-I) <= 1e-10
    end
    GB.CartesianShellification._collinear_terminal_geometry(direct,shell,slab,grid_axes,[-8.,0.,8.],3,3)
    @test (sum(coverage)-gap_count[]-exterior[],gap_count[],exterior[]) == (1029,98,98)
    @test all(==(1),coverage)
    @test_throws ArgumentError GB.CartesianShellification._collinear_terminal_geometry(direct,shell,slab,grid_axes,[-8.,0.,8.],3,8)
    for (z,Z) in [([],[]),([0.],[1.,1.]),([0.,0.],[1.,1.]),([NaN],[1.]),([0.],[-1.])]
        @test_throws ArgumentError build(z,Z)
    end
    for kw in [(core_side=true,), (core_side=2,), (angular_reference_count=0,),
        (outer_face_count=false,), (tail_spacing=Inf,), (transverse_spacing=0.,)]
        @test_throws ArgumentError build([-2.4,0.,2.4],ones(3);kw...)
    end
    @test_throws ArgumentError build([-.2,0.,.2],ones(3))
    @test_throws ArgumentError build([-2.4,0.,2.4],ones(3);outer_face_count=99, padding_transverse=6.)
end
