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
