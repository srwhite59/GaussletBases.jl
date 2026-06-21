using GaussletBases
using JLD2
using LinearAlgebra
using Test

const H_LOWEST = -0.49855234726272035
const H2_LOWEST = -0.79460371733658908
const H2_SELF_COULOMB = 0.4569117646737212
const ATOL = 1.0e-10
const PROVENANCE_KEYS = (
    :provenance_version, :producer, :route, :q, :core_spacing,
    :reference_spacing, :tail_spacing, :parent_axis_family, :parent_axis_counts,
    :mapping_kind, :mapping_d, :radius, :xmax_parallel, :xmax_transverse,
    :atom_symbols, :nuclear_charges, :atom_locations, :nup, :ndn,
    :final_dimension,
)

const H_BASIS = (;
    q = 5,
    core_spacing = 0.5,
    radius = 4.0,
    reference_spacing = 1.0,
    d = 0.3,
)

const H2_BASIS = (;
    q = 5,
    core_spacing = 0.5,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
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

function check_provenance_keys(file)
    for key in PROVENANCE_KEYS
        @test haskey(file, "producer_provenance/$(key)")
    end
end

@testset "public Cartesian base Hamiltonian" begin
    h = cartesian_base_hamiltonian(h_system(); basis = H_BASIS)
    h_lowest, _, _ = lowest_one_body(h)
    @test h isa CartesianIDAHamiltonian{Float64}
    @test h_lowest ≈ H_LOWEST atol = ATOL
    check_finite_symmetric(h.kinetic)
    check_finite_symmetric(h.electron_electron_ida)
    foreach(check_finite_symmetric, h.nuclear_attraction_unit_by_center)

    h2 = cartesian_base_hamiltonian(h2_system(); basis = H2_BASIS)
    h2_lowest, orbital, h2_one_body = lowest_one_body(h2)
    @test h2 isa CartesianIDAHamiltonian{Float64}
    @test size(h2.kinetic) == (471, 471)
    @test h2_lowest ≈ H2_LOWEST atol = ATOL
    @test self_coulomb(h2.electron_electron_ida, orbital) ≈ H2_SELF_COULOMB atol = ATOL
    check_finite_symmetric(h2.kinetic)
    check_finite_symmetric(h2.electron_electron_ida)
    foreach(check_finite_symmetric, h2.nuclear_attraction_unit_by_center)

    mktempdir() do dir
        h_path = joinpath(dir, "h_cartesian_ida.jld2")
        written_h = cartesian_base_hamiltonian(h_system(); basis = H_BASIS, hamfile = h_path)
        jldopen(h_path, "r") do file
            check_provenance_keys(file)
            @test file["producer_provenance/reference_spacing"] == 1.0
            @test file["producer_provenance/mapping_d"] == 0.3
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
            @test file["producer_provenance/mapping_kind"] === :multicenter_pqs_mapping
            @test file["producer_provenance/mapping_d"] === nothing
            @test file["producer_provenance/route"] === :z_axis_diatomic_pqs_base
        end
    end

    @test_throws ArgumentError cartesian_base_hamiltonian(
        merge(h_system(), (; extra = true)); basis = H_BASIS)
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h_system(); basis = merge(H_BASIS, (; xmax_parallel = 6.0)))
    @test_throws ArgumentError cartesian_base_hamiltonian(
        h_system(); basis = (; q = 5, core_spacing = 0.5, radius = 4.0))
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
end
