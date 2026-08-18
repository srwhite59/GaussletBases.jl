using GaussletBases
using LinearAlgebra

const SYSTEM = (;
    atom_symbols = ["H", "H"],
    nuclear_charges = [1.0, 1.0],
    atom_locations = [(0.0, 0.0, -1.0), (0.0, 0.0, 1.0)],
    nup = 1,
    ndn = 0,
)

const BASIS = (;
    ns = 4,
    core_spacing = 0.6,
    xmax_parallel = 3.0,
    xmax_transverse = 2.0,
    tail_spacing = 2.8,
)

const TOLERANCE = 1.0e-10

function build_h2plus(method)
    ham = cartesian_base_hamiltonian(
        SYSTEM;
        basis = merge(BASIS, (; nesting = method)),
    )
    h1 = one_body_hamiltonian(ham)
    vee = ham.electron_electron_ida
    eig = eigen(Symmetric(h1))
    energy = eig.values[1]
    orbital = eig.vectors[:, 1]
    residual = norm(h1 * orbital - energy * orbital)

    @assert ham isa CartesianIDAHamiltonian{Float64}
    @assert size(h1) == size(vee) == (293, 293)
    @assert ham.nup == 1 && ham.ndn == 0
    @assert ham.nuclear_charges == SYSTEM.nuclear_charges
    @assert ham.nuclear_positions == [0.0 0.0 -1.0; 0.0 0.0 1.0]
    @assert nuclear_repulsion(ham) == 0.5
    @assert all(isfinite, h1) && all(isfinite, vee)
    @assert norm(h1 - transpose(h1), Inf) <= TOLERANCE
    @assert norm(vee - transpose(vee), Inf) <= TOLERANCE
    @assert residual <= TOLERANCE
    return (; method, energy, residual)
end

for method in (:pqs, :wl)
    result = build_h2plus(method)
    println(
        result.method,
        ": dimension=293, lowest H1 energy=", result.energy,
        ", eigen-residual=", result.residual,
    )
end
