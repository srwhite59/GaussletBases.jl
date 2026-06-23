#!/usr/bin/env julia
# Canonical human-facing Cartesian producer. Edit these inputs directly, or pass
# one trusted local input file and optional key=value overrides.
using GaussletBases
# Public inputs
mode = :base                  # :base or :supplemented
Natom = 2                     # 1 atom, or 2 for homonuclear z-axis diatomic
R = 4.0                       # full bond length in bohr when Natom == 2
Z = 1.0
atom = "H"                    # label only; Z is the charge authority
nup = 1
ndn = 1
q = 5
core_spacing = 0.5
reference_spacing = 1.0
tail_spacing = 10.0
parent_axis_family = :G10
radius = 4.0
d = 0.3
xmax_parallel = 6.0
xmax_transverse = 4.0
basisname = "cc-pVTZ"
lmax = 1
uncontracted = false
width_filtering = nothing
basisfile = nothing
hamfile = "cartesian_ida_hamiltonian.jld2"
readback = true
public_inputs = (
    :mode, :Natom, :R, :Z, :atom, :nup, :ndn, :q, :core_spacing,
    :reference_spacing, :tail_spacing, :parent_axis_family, :radius, :d,
    :xmax_parallel, :xmax_transverse, :basisname, :lmax, :uncontracted,
    :width_filtering, :basisfile, :hamfile, :readback,
)
allowed_inputs = Set(public_inputs)
vars = Dict(name => getfield(Main, name) for name in public_inputs)
function apply_inputs!(vars, values, allowed_inputs)
    for (key, value) in pairs(values)
        name = Symbol(key)
        name in allowed_inputs || throw(ArgumentError("unknown driver input: $(name)"))
        vars[name] = value
    end
    return vars
end
args = collect(ARGS)
if !isempty(args) && !occursin("=", first(args))
    file_value = Base.include(Main, popfirst!(args))
    if file_value isa NamedTuple || file_value isa AbstractDict
        apply_inputs!(vars, file_value, allowed_inputs)
    else
        vars = Dict(name => getfield(Main, name) for name in public_inputs)
    end
end
for arg in args
    key, value = split(arg, "="; limit = 2)
    Symbol(key) in allowed_inputs || throw(ArgumentError("unknown driver input: $(key)"))
    vars[Symbol(key)] = Core.eval(Main, Meta.parse(value))
end
N = Int(vars[:Natom])
charge = Float64(vars[:Z])
label = String(vars[:atom])
system = if N == 1
    (;  atom_symbols = [label], nuclear_charges = [charge],
        atom_locations = [(0.0, 0.0, 0.0)],
        nup = vars[:nup], ndn = vars[:ndn])
elseif N == 2
    half_R = Float64(vars[:R]) / 2
    (;  atom_symbols = [label, label], nuclear_charges = [charge, charge],
        atom_locations = [(0.0, 0.0, -half_R), (0.0, 0.0, half_R)],
        nup = vars[:nup], ndn = vars[:ndn])
else
    throw(ArgumentError("Natom must be 1 or 2"))
end
common_basis = (; q = vars[:q], core_spacing = vars[:core_spacing],
    reference_spacing = vars[:reference_spacing],
    tail_spacing = vars[:tail_spacing],
    parent_axis_family = vars[:parent_axis_family])
basis = N == 1 ?
    (; common_basis..., radius = vars[:radius], d = vars[:d]) :
    (; common_basis...,
        xmax_parallel = vars[:xmax_parallel],
        xmax_transverse = vars[:xmax_transverse])
hamfile_value = String(vars[:hamfile])
isempty(hamfile_value) && throw(ArgumentError("hamfile must not be empty"))
elapsed = @elapsed ham = if vars[:mode] === :base
    GaussletBases.cartesian_base_hamiltonian(system; basis, hamfile = hamfile_value)
elseif vars[:mode] === :supplemented
    N == 2 || throw(ArgumentError("supplemented mode supports diatomics only"))
    supplement = (; basis_by_center = fill(String(vars[:basisname]), N),
        lmax = vars[:lmax], uncontracted = vars[:uncontracted],
        width_filtering = vars[:width_filtering], basisfile = vars[:basisfile])
    GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        system; basis, supplement, hamfile = hamfile_value)
else
    throw(ArgumentError("mode must be :base or :supplemented"))
end
dimension = size(ham.kinetic, 1)
if vars[:readback]
    loaded = GaussletBases.read_cartesian_ida_hamiltonian(hamfile_value)
    size(loaded.kinetic, 1) == dimension ||
        throw(ArgumentError("readback dimension mismatch"))
end
println("mode: ", vars[:mode])
println("Natom: ", N, "  atom: ", label, "  Z: ", charge)
println("electrons: nup=", vars[:nup], " ndn=", vars[:ndn])
println("hamfile: ", hamfile_value)
println("dimension: ", dimension)
println("elapsed_s: ", round(elapsed; digits = 3))
