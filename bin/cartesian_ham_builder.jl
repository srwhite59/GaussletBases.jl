#!/usr/bin/env julia
# Canonical human-facing Cartesian producer. Edit these inputs directly, or pass
# one trusted local input file and optional key=value overrides.
using GaussletBases

# Public inputs
gausslet_family = :G10         # parent gausslet family
Natom = 2                      # 1 atom, or 2 for homonuclear z-axis diatomic
atom = "H"                     # label only; Z is the charge authority
Z = 1.0
R = 4.0                        # full bond length in bohr when Natom == 2
nup = 1
ndn = 1
ns = 5                         # requested cube/source/nesting size
nesting = :pqs                  # :pqs or :wl construction family
source_span = :ordinary         # :ordinary or :mapped_comx for PQS source spans
core_spacing = 0.3             # near-nucleus spacing / atom mapping width
padding = 10.0                 # extra box padding beyond each nucleus, in bohr

basisname = nothing            # set e.g. "cc-pVTZ" for supplemented diatomic
lmax = 1                       # supplement angular cutoff
uncontracted = false           # false keeps contracted supplement functions
supplement_width_max = nothing # optional maximum supplement width
basisfile = nothing            # optional trusted local BasisSets path

hamfile = "cartesian_ida_hamiltonian.jld2"

# Run hooks
# Agent note: keep workflow controls here. If a Codex-run validation needs a new
# public run-level switch, ask for approval and add it here rather than bypassing
# this driver or adding route/internal controls below.
check_file = true              # read back output file to check it
print_contract = true
print_timing = true
expected_dimension = nothing

public_inputs = (
    :Natom, :R, :Z, :atom, :nup, :ndn, :ns, :core_spacing, :padding,
    :gausslet_family, :nesting, :source_span, :basisname, :lmax, :uncontracted,
    :supplement_width_max, :basisfile, :hamfile, :check_file,
    :print_contract, :print_timing, :expected_dimension,
)

# Apply trusted inputs
vars = Dict(name => getfield(Main, name) for name in public_inputs)
function apply_inputs!(values)
    for (key, value) in pairs(values)
        name = Symbol(key)
        name in public_inputs || throw(ArgumentError("unknown driver input: $(name)"))
        vars[name] = value
    end
end
args = collect(ARGS)
if !isempty(args) && !occursin("=", first(args))
    input_path = abspath(popfirst!(args))
    file_value = Base.include(Main, input_path)
    file_value isa NamedTuple || file_value isa AbstractDict ?
        apply_inputs!(file_value) :
        (vars = Dict(name => getfield(Main, name) for name in public_inputs))
end
for arg in args
    key, value = split(arg, "="; limit = 2)
    Symbol(key) in public_inputs || throw(ArgumentError("unknown driver input: $(key)"))
    vars[Symbol(key)] = Core.eval(Main, Meta.parse(value))
end

function print_terminal_inventory(base)
    inventory = hasproperty(base, :terminal_inventory) ? base.terminal_inventory : nothing
    isnothing(inventory) && return nothing
    println("terminal inventory: base_final_dimension=", inventory.final_dimension)
    println("  region kind lowering support final ratio class slab")
    for row in inventory.rows
        ratio = round(row.compression_ratio; digits = 3)
        slab = row.slab_axis === :unavailable ? "" :
            " slab=$(row.slab_axis)/$(row.slab_side)/t=$(row.slab_thickness)/$(row.slab_stack_index)-$(row.slab_stack_count)"
        println("  ", row.region_key, " ", row.region_kind, " ", row.lowering_kind, " support=", row.support_rows, " final=", row.final_cols, " ratio=", ratio, " ", row.realization_class, slab)
    end
    return nothing
end

total_start = time()
# Construct public contract
contract_start = time()
N = Int(vars[:Natom])
charge = Float64(vars[:Z])
label = String(vars[:atom])
nesting_value = Symbol(vars[:nesting])
nesting_value in (:pqs, :wl) || throw(ArgumentError("nesting must be :pqs or :wl"))
source_span_value = Symbol(vars[:source_span])
source_span_value in (:ordinary, :mapped_comx) || throw(ArgumentError("source_span must be :ordinary or :mapped_comx"))
N in (1, 2) || throw(ArgumentError("Natom must be 1 or 2"))
system = if N == 1
    (;  atom_symbols = [label], nuclear_charges = [charge], atom_locations = [(0.0, 0.0, 0.0)], nup = vars[:nup], ndn = vars[:ndn])
else
    half_R = Float64(vars[:R]) / 2
    (;  atom_symbols = [label, label], nuclear_charges = [charge, charge], atom_locations = [(0.0, 0.0, -half_R), (0.0, 0.0, half_R)], nup = vars[:nup], ndn = vars[:ndn])
end
common_basis = (; ns = vars[:ns], core_spacing = vars[:core_spacing],
    parent_axis_family = vars[:gausslet_family], nesting = nesting_value,
    source_span = source_span_value)
basis = N == 1 ? (; common_basis..., radius = vars[:padding]) :
    (; common_basis..., xmax_parallel = Float64(vars[:R]) / 2 + vars[:padding],
        xmax_transverse = vars[:padding])
hamfile_value = String(vars[:hamfile])
isempty(hamfile_value) && throw(ArgumentError("hamfile must not be empty"))
supplemented = !isnothing(vars[:basisname])
supplement = nothing
if supplemented
    width_filtering = isnothing(vars[:supplement_width_max]) ? nothing : (; max_width = vars[:supplement_width_max])
    supplement = (; basis_by_center = fill(String(vars[:basisname]), N), lmax = vars[:lmax], uncontracted = vars[:uncontracted], width_filtering, basisfile = vars[:basisfile])
end
contract_elapsed = time() - contract_start

# Review public contract
if vars[:print_contract]
    println("system: Natom=", N, " atom=", label, " Z=", charge, " R=", vars[:R], " nup=", vars[:nup], " ndn=", vars[:ndn])
    println("basis: ns=", vars[:ns], " nesting=", nesting_value, " source_span=", source_span_value, " core_spacing=", vars[:core_spacing], " padding=", vars[:padding], " gausslet_family=", vars[:gausslet_family])
    !isnothing(supplement) && println("supplement: basisname=", vars[:basisname], " lmax=", vars[:lmax], " uncontracted=", vars[:uncontracted], " supplement_width_max=", vars[:supplement_width_max], " basisfile=", vars[:basisfile])
    println("hamfile: ", hamfile_value)
    println("hooks: check_file=", vars[:check_file], " print_contract=", vars[:print_contract], " print_timing=", vars[:print_timing], " expected_dimension=", vars[:expected_dimension])
end

# Build Hamiltonian
stage_timings = Pair{String,Float64}[]
build_start = time()

stage_start = time()
base = GaussletBases.cartesian_base_working_basis(system; basis, supplemented)
push!(stage_timings, "base working basis" => time() - stage_start)
vars[:print_contract] && print_terminal_inventory(base)

stage_start = time()
base_products = GaussletBases.cartesian_base_products(base)
push!(stage_timings, "base products" => time() - stage_start)

stage_start = time()
base_unit_nuclear = GaussletBases.cartesian_base_unit_nuclear(base)
push!(stage_timings, "base unit nuclear" => time() - stage_start)

stage_start = time()
base_vee = GaussletBases.cartesian_base_vee(base)
push!(stage_timings, "base electron-electron" => time() - stage_start)

stage_start = time()
base_ham = GaussletBases.cartesian_base_hamiltonian_assembly(
    base, base_products, base_unit_nuclear, base_vee;
    hamfile = supplemented ? nothing : hamfile_value)
push!(stage_timings, "base Hamiltonian" => time() - stage_start)

if supplemented
    stage_start = time()
    supplement_basis = GaussletBases.cartesian_residual_gto_supplement_basis(base, supplement)
    push!(stage_timings, "supplement basis" => time() - stage_start)

    stage_start = time()
    residual = GaussletBases.cartesian_residual_gto_augmentation(base, supplement_basis)
    push!(stage_timings, "residual augmentation" => time() - stage_start)

    stage_start = time()
    augmented_products = GaussletBases.cartesian_residual_gto_augmented_products(
        base, supplement_basis, residual; base_kinetic = base_ham.kinetic)
    push!(stage_timings, "augmented products" => time() - stage_start)

    stage_start = time()
    augmented_unit_nuclear = GaussletBases.cartesian_residual_gto_augmented_unit_nuclear(
        base, residual, augmented_products;
        base_unit_nuclear = base_ham.nuclear_attraction_unit_by_center)
    push!(stage_timings, "augmented unit nuclear" => time() - stage_start)

    stage_start = time()
    augmented_vee = GaussletBases.cartesian_residual_gto_augmented_vee(
        base, base_ham, residual, augmented_products, augmented_unit_nuclear)
    push!(stage_timings, "augmented electron-electron" => time() - stage_start)

    stage_start = time()
    ham = GaussletBases.cartesian_residual_gto_mwg_hamiltonian_assembly(
        base, base_ham, supplement_basis, residual, augmented_products,
        augmented_unit_nuclear, augmented_vee;
        hamfile = hamfile_value)
    push!(stage_timings, "supplemented Hamiltonian" => time() - stage_start)
else
    ham = base_ham
end
build_elapsed = time() - build_start
dimension = size(ham.kinetic, 1)

# Check artifact
check_start = time()
isnothing(vars[:expected_dimension]) || dimension == Int(vars[:expected_dimension]) ||
    throw(ArgumentError("dimension mismatch"))
if vars[:check_file]
    loaded = GaussletBases.read_cartesian_ida_hamiltonian(hamfile_value)
    size(loaded.kinetic, 1) == dimension ||
        throw(ArgumentError("file check dimension mismatch"))
end
check_elapsed = time() - check_start

# Summarize run
total_elapsed = time() - total_start
println("supplemented: ", supplemented)
println("dimension: ", dimension)
if vars[:print_timing]
    println("timing_s:")
    println("  contract=", round(contract_elapsed; digits = 3))
    for (label, elapsed) in stage_timings
        println("  ", label, "=", round(elapsed; digits = 3))
    end
    println("  build total=", round(build_elapsed; digits = 3))
    println("  check=", round(check_elapsed; digits = 3))
    println("  total=", round(total_elapsed; digits = 3))
end
